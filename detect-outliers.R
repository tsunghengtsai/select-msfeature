# Load required packages --------------------------------------------------
library(tidyverse)
library(stringr)
library(broom)


# Load functions and dataset ----------------------------------------------
source("utilities.R")

# DDA datasets 
load("data/DDA/iPRG_20170418/iprg.mq.RData")
df_mss <- iprg.mq

# load("data/DDA/iPRG_20170418/iprg.progenesis.RData")
# df_mss <- iprg.progenesis

# load("data/DDA/iPRG_20170418/iprg.skyline.RData")
# df_mss <- iprg.skyline

# DIA datasets


# Data wrangling ----------------------------------------------------------
# Convert to the working format
df_allftr <- pdata2ftr(df_mss)


# For each protein, remove uncovered runs and undetectable features 
# (presumed by subsequent analysis)
df_allftr <- df_allftr %>% 
    group_by(protein, run) %>% filter(any(is_obs)) %>% ungroup() %>% 
    group_by(protein, feature) %>% filter(any(is_obs)) %>% ungroup()


# Residuals around TMP estimates ------------------------------------------
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

# Technical variability
nested_prot <- nested_prot %>% 
    mutate(
        mp_fit = map(subdata, medpolish_df), 
        mp_resid = map(mp_fit, residuals_mp), 
        mad_res = map_dbl(mp_resid, ~ mad(.$log2res, na.rm = T)), 
        # nb_obs = map_int(mp_resid, ~ sum(!is.na(.$log2res))), 
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


# Save .RDS files ---------------------------------------------------------
# saveRDS(as.data.frame(df_rmftr), "output/iprg.mq.rmftr.RDS")
# saveRDS(as.data.frame(df_rmftr), "output/iprg.progenesis.rmftr.RDS")
# saveRDS(as.data.frame(df_rmftr), "output/iprg.skyline.rmftr.RDS")

# saveRDS(nested_prot, "output/iprg.mq.nest.RDS")
# saveRDS(nested_prot, "output/iprg.progenesis.nest.RDS")
# saveRDS(nested_prot, "output/iprg.skyline.nest.RDS")


# Print profiles to file --------------------------------------------------
# Features to be removed
# pdf("output/iprg_mq_rm_04.pdf", width = 9, height = 6)
# pdf("output/iprg_progenesis_rm.pdf", width = 9, height = 6)
pdf("output/iprg_skyline_rm01.pdf", width = 9, height = 6)
for (i in which(nested_prot$protein %in% df_rmftr$protein)[1:500]) {
    print(plot_profile_nest(nested_prot, i))
}
dev.off()


# Spike-in
pdf("output/iprg_mq_spike.pdf", width = 9, height = 6)
pdf("output/iprg_progenesis_spike.pdf", width = 9, height = 6)
pdf("output/iprg_skyline_spike.pdf", width = 9, height = 6)
for (i in which(str_detect(nested_prot$protein, "P44015|P55752|P44374|P44983|P44683|P55249"))) {
    print(plot_profile_nest(nested_prot, i))
}
dev.off()


# pdf("output/iprg_progenesis_02.pdf", width = 9, height = 6)
# pdf("output/iprg_skyline_02.pdf", width = 9, height = 6)
# for (i in 101:200) {
#     print(plot_profile_nest(nested_prot, i))
# }
# dev.off()

# Initial attempt for shrinking variance estimates ------------------------
# fit_gamma <- function(x) {
#     mu_x <- mean(x)
#     var_x <- var(x)
#     shape_g <- mu_x ^ 2 / var_x
#     rate_g <- mu_x / var_x
#     fit_g <- MASS::fitdistr(x, densfun = "gamma", 
#                             start = list(shape = shape_g, rate = rate_g))
#     
#     return(fit_g$estimate)
# }
# 
# prot_res <- nested_prot %>% 
#     select(protein, mad_res, nb_obs) %>% 
#     mutate(mad_res2 = mad_res ^ 2)
# 
# mad2_sample <- prot_res %>% filter(mad_res != 0, nb_obs > 60) %>% .$mad_res2 
# param_gamma <- fit_gamma(mad2_sample)
# hist(mad2_sample, pch=20, breaks=100, prob=TRUE, main="")
# curve(dgamma(x, param_gamma[1], param_gamma[2]), col="red", lwd=2, add=T)

