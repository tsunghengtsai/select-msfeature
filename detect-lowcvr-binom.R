# Load required packages --------------------------------------------------
library(MASS)
library(tidyverse)
library(stringr)
library(broom)
library(VGAM)


# Load functions and dataset ----------------------------------------------

# DDA datasets 
# load("data/DDA_iPRG/dda.iprg.mq.rda")
# load("data/DDA_iPRG/dda.iprg.pg.rda")
# load("data/DDA_iPRG/dda.iprg.sl.rda")
# 
# load("data/DDA_Cox/dda.cox.mq.rda")
# load("data/DDA_Cox/dda.cox.pg.rda")
# load("data/DDA_Cox/dda.cox.sl.rda")
# 
# load("data/DDA_Choi/dda.choi.mq.rda")
# load("data/DDA_Choi/dda.choi.pg.rda")
# load("data/DDA_Choi/dda.choi.pd.rda")
load("data/DDA_Choi/dda.choi.sl.rda")

# DIA datasets
# load("data/DIA_Bruderer/dia.bruderer.sl.rda")
# load("data/DIA_Bruderer/dia.bruderer.sn.rda")

# load("data/DIA_Navarro/dia.navarro.du.rda")
# load("data/DIA_Navarro/dia.navarro.os.rda")
# load("data/DIA_Navarro/dia.navarro.sl.rda")
# load("data/DIA_Navarro/dia.navarro.sn.rda")

# load("data/DIA_Eye/dia.eye.sn.rda")
# load("data/DIA_Selevsek/dia.selevsek.sn.rda")

# SRM datasets



# Housekeeping ------------------------------------------------------------
dset_name <- ls()
df_mss <- eval(parse(text = dset_name))

source("utilities.R")


# Data wrangling 
# Convert to the working format
df_allftr <- pdata2ftr(df_mss)

# For each protein, remove uncovered runs and undetectable features 
# (presumed by subsequent analysis)
df_allftr <- df_allftr %>% 
    group_by(protein, run) %>% filter(any(is_obs)) %>% ungroup() %>% 
    group_by(protein, feature) %>% filter(any(is_obs)) %>% ungroup()


# Proportion of observation -----------------------------------------------

obs_protein <- df_allftr %>% 
    group_by(protein) %>% 
    summarise(
        nb_run = n_distinct(run), 
        nb_feature = n_distinct(feature), 
        nb_obs = sum(is_obs)
    ) %>% 
    mutate(nb_full = nb_run * nb_feature, p_obs = nb_obs / nb_full) %>% 
    mutate(min_obs = qbinom(0.05, nb_run, p_obs))

nested_prot <- df_allftr %>% 
    select(protein, run, peptide, feature, log2inty, is_obs) %>% 
    group_by(protein) %>% 
    nest()

nested_prot <- nested_prot %>% 
    left_join(obs_protein %>% select(protein, min_obs))


# Coverage
list_data <- vector("list", length = nrow(nested_prot))
for (i in seq_along(list_data)) {
    list_data[[i]] <- flag_lowcover(nested_prot$data[[i]]) %>% 
        group_by(feature) %>% 
        mutate(is_lowbinom = sum(is_obs) < nested_prot$min_obs[i]) %>% 
        ungroup()
}
nested_prot$data <- list_data

nested_prot <- nested_prot %>% 
    mutate(cover = map(data, ~distinct(., feature, is_lowcvr, is_lowbinom)))

nested_ftr <- nested_prot %>% 
    select(protein, cover) %>% 
    unnest() %>% 
    filter(is_lowbinom) %>% 
    nest(-protein)


# Profile -----------------------------------------------------------------

plot_profile_ftr <- function(nested, idx, ftr, yrange = c(10, 35)) {
    oneprot <- nested$data[[idx]] %>%
        mutate(hivar = if_else(feature %in% ftr, "high var / low coverage", "low var"))
    oneprot %>%
        ggplot(aes(run, log2inty, color = peptide, group = feature)) +
        geom_point(size = 3, alpha = 0.5) + geom_line() +
        scale_alpha_discrete(range = c(0.5, 1)) +
        facet_wrap(~ hivar, drop = F) +
        ggtitle(nested$protein[idx]) +
        coord_cartesian(ylim = yrange) +
        theme(legend.position = "none", axis.text.x = element_blank())
}

pdf(paste0(str_replace_all(dset_name, "\\.", "_"), "_binom.pdf"), width = 9, height = 6)
for (i in 1:min(nrow(nested_ftr), 200)) {
    idx <- which(nested_prot$protein == nested_ftr$protein[i])
    print(plot_profile_ftr(nested_prot, idx, nested_ftr$data[[i]]$feature, yrange = c(10, 35)))
}
dev.off()


# spike <- "P44015|P55752|P44374|P44983|P44683|P55249"
spike <- "SS1_P00432|SS1_P00711|SS1_P02701|SS1_P02754|SS1_P0CG53|SS1_P24627|SS1_P68082|SS1_P80025|SS1_Q29443|SS2_P00915|SS2_P00921|SS2_P01008|SS2_P01012|SS2_P02662|SS2_P02663|SS2_P02666|SS2_P02787|SS2_P05307|SS2_P61769|SS3_P00563|SS3_P00698|SS3_P02769|SS3_Q3SX14|SS3_Q58D62|SS4_P00004|SS4_P00442|SS4_P01133|SS4_P02753"
# spike <- "P02754|P00921|P80025|P02662|P00366|P12799|P02672|P02789|P02676|P61823|P68082|P02666"
pdf(paste0(str_replace_all(dset_name, "\\.", "_"), "_spike_binom.pdf"), width = 9, height = 6)
for (i in which(str_detect(nested_prot$protein, spike))) {
    if (nested_prot$protein[i] %in% nested_ftr$protein) {
        ii <- which(nested_ftr$protein == nested_prot$protein[i])
        print(plot_profile_ftr(nested_prot, i, nested_ftr$data[[ii]]$feature, yrange = c(10, 35)))
    }
}
dev.off()





# EB shrinkage toward a beta ----------------------------------------------

tmp <- df_allftr %>% 
    group_by(protein) %>% 
    summarise(
        nb_run = n_distinct(run), 
        nb_feature = n_distinct(feature), 
        nb_obs = sum(is_obs)
    ) %>% 
    mutate(nb_full = nb_run * nb_feature, p_obs = nb_obs / nb_full)

samp_mean <- mean(tmp$p_obs)
samp_var <- var(tmp$p_obs)
a <- samp_mean * (samp_mean * (1 - samp_mean) / samp_var - 1)
b <- (1 - samp_mean) * (samp_mean * (1 - samp_mean) / samp_var - 1)

# m <- MASS::fitdistr(tmp$p_obs[tmp$nb_feature > 5 & tmp$p_obs < 0.99], dbeta, start = list(shape1=a, shape2=b))

# negative log likelihood of data given alpha; beta
ll <- function(alpha, beta) {
    -sum(dbetabinom.ab(tmp$nb_obs, tmp$nb_full, alpha, beta, log = TRUE))
}
m <- mle(ll, start = list(alpha = a, beta = b), method = "L-BFGS-B")
coef(m)


# Plots -------------------------------------------------------------------

tmp %>% 
    # filter(nb_feature > 5) %>% 
    ggplot(aes(p_obs)) + 
    geom_histogram(aes(y = ..density..), binwidth = 0.02) + 
    stat_function(fun = dbeta, args = list(shape1 = coef(m)[1], shape2 = coef(m)[2]), col = 'red')


# tmp %>%
#     filter(nb_feature > 5, p_obs < 0.99) %>%
#     ggplot(aes(p_obs)) +
#     geom_histogram(aes(y = ..density..), binwidth = 0.02) +
#     stat_function(fun = dbeta, args = list(shape1 = m$estimate[1], shape2 = m$estimate[2]), col = 'red')

tmp <- tmp %>% 
    mutate(p_adj = (nb_obs + coef(m)[1]) / (nb_full + coef(m)[1] + coef(m)[2]))
    # mutate(p_adj = (nb_obs + m$estimate[1]) / (nb_full + m$estimate[1] + m$estimate[2]))

tmp %>% 
    ggplot(aes(p_obs, p_adj, color = nb_full)) + 
    geom_point() + 
    geom_abline(color = "red") + 
    geom_hline(yintercept = coef(m)[1] / (coef(m)[1] + coef(m)[2]), color = "red", linetype = "dashed")
    # geom_hline(yintercept = m$estimate[1] / (m$estimate[1] + m$estimate[2]), color = "red", linetype = "dashed")


binom.test(4, 15, 0.6, alternative = "less")

