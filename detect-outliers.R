# Load required packages --------------------------------------------------
library(MASS)
library(tidyverse)
library(stringr)
library(broom)


# Load functions and dataset ----------------------------------------------

# DDA datasets 
load("data/DDA/iPRG_20170418/iprg.mq.RData")
# load("data/DDA/iPRG_20170418/iprg.progenesis.RData")
# load("data/DDA/iPRG_20170418/iprg.skyline.RData")

# DIA datasets
# load("data/DIA/Bruderer_20170418/bruderer.skyline.RData")
# load("data/DIA/Bruderer_20170418/bruderer.SN.RData")

dset_name <- ls()
df_mss <- eval(parse(text = dset_name))

source("utilities.R")


# Data wrangling ----------------------------------------------------------
# Convert to the working format
df_allftr <- pdata2ftr(df_mss)


# For each protein, remove uncovered runs and undetectable features 
# (presumed by subsequent analysis)
df_allftr <- df_allftr %>% 
    group_by(protein, run) %>% filter(any(is_obs)) %>% ungroup() %>% 
    group_by(protein, feature) %>% filter(any(is_obs)) %>% ungroup()


# Residuals around RLM estimates ------------------------------------------
nested_prot <- df_allftr %>% 
    select(protein, run, peptide, feature, log2inty, is_obs) %>% 
    group_by(protein) %>% 
    nest()

# Coverage
nested_prot <- nested_prot %>% 
    mutate(
        data = map(data, flag_lowcover), 
        subdata = map(data, ~ filter(., !is_lowcvr))
    )

# Robust linear model
nested_prot <- nested_prot %>% 
    mutate(
        rlm_fit = map(
            subdata, 
            possibly(~ rlm(log2inty ~ run + feature, data = ., scale.est = "Huber", 
                           na.action = na.exclude), NULL)), 
        dfree = map_dbl(subdata, approx_df2)
    ) %>% 
    mutate(s_res = map_dbl(rlm_fit, ~ ifelse(!is_null(.), summary(.)$sigma, NA)))


# Shrinkage variance estimation with limma --------------------------------

prot_res <- nested_prot %>% 
    select(protein, s_res, dfree) %>% 
    mutate(var_res = s_res ^ 2)

prot_res2 <- prot_res %>% filter(!is.na(s_res), dfree > 0)
eb_fit <- limma::squeezeVar(var = prot_res2$var_res, df = prot_res2$dfree)
ebrob_fit <- limma::squeezeVar(var = prot_res2$var_res, df = prot_res2$dfree, robust = TRUE)

prot_res2 <- prot_res2 %>% 
    mutate(var_eb = eb_fit$var.post, var_ebrob = ebrob_fit$var.post)


# variance estimates before/after shrinkage -------------------------------
prot_res2 %>% ggplot(aes(var_res, var_eb, colour = dfree)) + 
    geom_point() + geom_abline(color = "red")

# EB fit
prot_res2 %>% ggplot(aes(sqrt(var_res), sqrt(var_eb), colour = dfree)) + 
    geom_point() + 
    geom_abline(color = "red") + 
    geom_hline(yintercept = sqrt(eb_fit$var.prior), color = "red", lty = 2) +
    scale_colour_gradient(trans = "log", breaks = 4 ^ (1:5)) + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

# Robust EB fit
prot_res2 %>% ggplot(aes(sqrt(var_res), sqrt(var_ebrob), colour = dfree)) + 
    geom_point() + 
    geom_abline(color = "red") + 
    geom_hline(yintercept = sqrt(ebrob_fit$var.prior), color = "red", lty = 2) +
    scale_colour_gradient(trans = "log", breaks = 4 ^ (1:5)) + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

# Difference between EB fit and original scale estimate
prot_res2 %>% mutate(d_scale = sqrt(var_ebrob) - sqrt(var_res)) %>% 
    ggplot(aes(sqrt(var_res), d_scale, colour = dfree)) + 
    geom_point() + 
    geom_hline(yintercept = 0) + 
    scale_colour_gradient(trans = "log", breaks = 4 ^ (1:5))

prot_res2 %>% mutate(d_mad = sqrt(var_ebrob) - sqrt(var_res)) %>% 
    filter(d_mad > 0.1, dfree > 5)

# Histogram of all variance estimates
prot_res2 %>% ggplot(aes(var_res)) + 
    geom_histogram(binwidth = 0.01) + 
    geom_vline(xintercept = (ebrob_fit$var.prior), color = "red", lty = 2) + 
    coord_cartesian(xlim = c(0, 0.5))

# Compare histogram of variances (with the same df) & fitted model
fd_fit <- limma::fitFDistRobustly(x = prot_res2$var_res, df1 = prot_res2$dfree)
sdf <- function(x, df1, df2, s) {df(x / s, df1, df2) / s}  # scaled F distribution

prot_res2 %>% filter(dfree == 11) %>% 
    ggplot(aes(var_res)) + 
    geom_histogram(aes(y = ..density..), binwidth = 0.025) + 
    stat_function(fun = sdf, args = list(df1 = 11, df2 = fd_fit$df2, s = fd_fit$scale), col = 'red')


# Flag outliers with scale estimates --------------------------------------

nested_prot2 <- nested_prot %>% 
    inner_join(prot_res2) %>% 
    mutate(subdata = map2(subdata, rlm_fit, ~ mutate(.x, log2res = residuals(.y))))

# Flag outliers based on robust EB estimates (other options: EB, MAD)
nested_prot2 <- nested_prot2 %>% 
    mutate(subdata = map2(subdata, var_ebrob, ~ mutate(.x, is_olr = abs(log2res) > 3 * sqrt(.y))))

# Flag highly-variable features
nested_prot2 <- nested_prot2 %>% 
    mutate(subdata = map(subdata, flag_hivar))

# Features to be removed
nested_prot2 <- nested_prot2 %>% 
    mutate(feature_rm = map2(data, subdata, list_rmfeature))

df_rmftr <- nested_prot2 %>% select(protein, feature_rm) %>% unnest(feature_rm)

# Peaks (feature in run) to be removed
df_rmpk <- nested_prot2 %>% 
    mutate(outlier = map(subdata, ~ select(filter(., is_olr), run, feature))) %>% 
    select(protein, outlier) %>% 
    unnest(outlier)

# df_allftr %>% anti_join(rename(df_rmftr, feature = feature_rm))
# df_allftr %>% anti_join(rename(df_rmftr, feature = feature_rm)) %>% anti_join(df_rmpk)


# Return features and peaks for removal -----------------------------------

# Save .RDS files
saveRDS(as.data.frame(df_rmftr), paste0("output/", dset_name, ".rmftr.RDS"))
saveRDS(as.data.frame(df_rmpk), paste0("output/", dset_name, ".rmpk.RDS"))


# Print profiles to file --------------------------------------------------
# Features to be removed
# pdf("output/iprg_mq_rm_01.pdf", width = 9, height = 6)
# pdf("output/iprg_progenesis_rm.pdf", width = 9, height = 6)
# pdf("output/iprg_skyline_rm01.pdf", width = 9, height = 6)
pdf("output/bruderer_skyline_rm01.pdf", width = 9, height = 6)
# pdf("output/bruderer_SN_rm01.pdf", width = 9, height = 6)
for (i in which(nested_prot2$protein %in% df_rmftr$protein)[1:500]) {
    if (!is.na(i)) {print(plot_profile_nest(nested_prot2, i, yrange = c(5, 25)))}
}
dev.off()


# Spike-in
# pdf("output/iprg_mq_spike.pdf", width = 9, height = 6)
# pdf("output/iprg_progenesis_spike.pdf", width = 9, height = 6)
pdf("output/iprg_skyline_spike.pdf", width = 9, height = 6)
for (i in which(str_detect(nested_prot2$protein, "P44015|P55752|P44374|P44983|P44683|P55249"))) {
    print(plot_profile_nest(nested_prot2, i))
}
dev.off()


# Median polish -----------------------------------------------------------

# Scale estimate through median polish is not well calibrated
nested_prot <- nested_prot %>%
    mutate(
        mp_fit = map(subdata, medpolish_df),
        mp_resid = map(mp_fit, residuals_mp),
        mad_res = map_dbl(mp_resid, ~ mad(.$log2res, na.rm = T)),
        dfree = map_dbl(subdata, approx_df2_mad)
    )


# Shrinkage variance estimation with limma --------------------------------
prot_res <- nested_prot %>% 
    select(protein, mad_res, dfree) %>% 
    mutate(var_res = mad_res ^ 2)

prot_res2 <- prot_res %>% filter(mad_res != 0, dfree > 0)
eb_fit <- limma::squeezeVar(var = prot_res2$var_res, df = prot_res2$dfree)
reb_fit <- limma::squeezeVar(var = prot_res2$var_res, df = prot_res2$dfree, robust = TRUE)

prot_res2 <- prot_res2 %>% 
    mutate(var_eb = eb_fit$var.post, var_reb = reb_fit$var.post)


# variance estimates before/after shrinkage -------------------------------
prot_res2 %>% ggplot(aes(var_res, var_eb, colour = dfree)) + 
    geom_point() + geom_abline(color = "red")

# EB fit
prot_res2 %>% ggplot(aes(sqrt(var_res), sqrt(var_eb), colour = dfree)) + 
    geom_point() + 
    geom_abline(color = "red") + 
    geom_hline(yintercept = sqrt(eb_fit$var.prior), color = "red", lty = 2) +
    scale_colour_gradient(trans = "log", breaks = 4 ^ (1:5)) + 
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6))

# Robust EB fit
prot_res2 %>% ggplot(aes(sqrt(var_res), sqrt(var_reb), colour = dfree)) + 
    geom_point() + 
    geom_abline(color = "red") + 
    geom_hline(yintercept = sqrt(reb_fit$var.prior), color = "red", lty = 2) +
    scale_colour_gradient(trans = "log", breaks = 4 ^ (1:5)) + 
    coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6))

prot_res2 %>% mutate(d_mad = sqrt(var_reb) - sqrt(var_res)) %>% 
    ggplot(aes(sqrt(var_res), d_mad, colour = dfree)) + 
    geom_point() + 
    geom_hline(yintercept = 0) + 
    scale_colour_gradient(trans = "log", breaks = 4 ^ (1:5))

prot_res2 %>% ggplot(aes(var_res)) + 
    geom_histogram(binwidth = 0.005) + 
    geom_vline(xintercept = (reb_fit$var.prior), color = "red", lty = 2)


# Flag outliers with scale estimates --------------------------------------
nested_prot <- nested_prot %>% left_join(prot_res2)

# Flag outliers based on robust EB estimates (other options: EB, MAD)
nested_prot <- nested_prot %>% 
    mutate(mp_resid = map2(mp_resid, var_reb, ~ mutate(.x, is_olr = abs(log2res) > 3 * sqrt(.y))))

# Flag highly-variable features
nested_prot <- nested_prot %>% 
    mutate(mp_resid = map(mp_resid, flag_hivar))

# Features to be removed
nested_prot <- nested_prot %>% 
    mutate(feature_rm = map2(data, mp_resid, list_rmfeature))

df_rmftr <- nested_prot %>% select(protein, feature_rm) %>% unnest(feature_rm)

