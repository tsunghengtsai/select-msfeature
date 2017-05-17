# Load required packages --------------------------------------------------
library(MASS)
library(tidyverse)
library(stringr)
library(broom)


# Load functions and dataset ----------------------------------------------

# DDA datasets 
# load("data/DDA/iPRG_20170418/iprg.mq.RData")
# load("data/DDA/iPRG_20170418/iprg.progenesis.RData")
# load("data/DDA/iPRG_20170418/iprg.skyline.RData")

# DIA datasets
# load("data/DIA/Bruderer_20170418/bruderer.skyline.RData")
load("data/DIA/Bruderer_20170418/bruderer.SN.RData")

dset_name <- ls()
df_mss <- eval(parse(text = dset_name))

source("utilities.R")


# Functions for feature variance ------------------------------------------

# Sample variance and size per feature
ftrsum <- function(fit) {
    fs <- augment(fit) %>% 
        group_by(feature) %>% 
        summarise(v = var(.resid), n = n())
    
    return(fs)
}

# Posterior probability of two-component mixture model
calc_mixprob <- function(x, mix_mu, mix_sigma, mix_lambda) {
    logodd <- dnorm(x, mix_mu[2], mix_sigma[2], log = TRUE) - 
        dnorm(x, mix_mu[1], mix_sigma[1], log = TRUE) + 
        log(mix_lambda[2]) - log(mix_lambda[1])
    
    return(1 / (1 + exp(logodd)))
}

plot_profile_ftr <- function(nested, idx, ftr, yrange = c(10, 35)) {
    oneprot <- nested$subdata[[idx]] %>% 
        mutate(hivar = if_else(feature %in% ftr, "high var", "low var"))
    oneprot %>% 
        ggplot(aes(run, log2inty, color = peptide, group = feature)) + 
        geom_point(size = 3, alpha = 0.5) + geom_line() + 
        scale_alpha_discrete(range = c(0.5, 1)) + 
        facet_wrap(~ hivar, drop = F) + 
        ggtitle(nested$protein[idx]) + 
        coord_cartesian(ylim = yrange) + 
        theme(legend.position = "none", axis.text.x = element_blank())
}


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
        rlm_fit = map(subdata, possibly(~ rlm(log2inty ~ run + feature, data = ., scale.est = "Huber"), NULL)), 
        dfree = map_dbl(subdata, approx_df2)
    ) %>% 
    mutate(s_res = map_dbl(rlm_fit, ~ ifelse(!is_null(.), summary(.)$sigma, NA)))


# High-var features -------------------------------------------------------

# Sample variance of feature
nested_prot2 <- nested_prot %>% 
    filter(dfree > 0) %>% 
    mutate(ftr = map(rlm_fit, ftrsum))

ftrvar <- nested_prot2 %>% unnest(ftr)

# Empirical Bayes
ebrob_fit <- limma::squeezeVar(var = ftrvar$v, df = ftrvar$n - 1, robust = TRUE)
eb_fit <- limma::squeezeVar(var = ftrvar$v, df = ftrvar$n - 1)
ftrvar <- ftrvar %>% 
    mutate(
        var_eb = eb_fit$var.post, 
        var_ebrob = ebrob_fit$var.post
    )

ftrvar <- ftrvar %>% 
    filter(v != 0) %>% 
    mutate(lv = log(v), lv_eb = log(var_eb), lv_ebrob = log(var_ebrob))

# ftrvar %>% ggplot(aes(lv)) + geom_histogram(aes(y = ..density..), binwidth = 0.025)
# ftrvar %>% ggplot(aes(lv_eb)) + geom_histogram(aes(y = ..density..), binwidth = 0.025)

# Mixture model
library(mixtools)
mix <- normalmixEM(ftrvar$lv_eb, mu = c(-3, 0))
# plot(mix, density = TRUE)
# summary(mix)

ftrvar <- ftrvar %>% 
    mutate(mixprob = calc_mixprob(lv_eb, mix$mu, mix$sigma, mix$lambda))

ftr_hivar <- ftrvar %>% filter(mixprob <= 0.05)
nested_ftr <- ftr_hivar %>% nest(-protein)


# Assuming most are from one single distribution
# ftrvar_sum <- ftrvar %>% 
#     filter(v != 0) %>% 
#     mutate(lv = log(v)) %>% 
#     summarise(lv_med = median(lv), lv_mad = mad(lv))
# nested_ftr <- ftr_hivar %>% nest(-protein)

# ftrvar %>% ggplot(aes(lv)) + 
#     geom_histogram(aes(y = ..density..), binwidth = 0.025) + 
#     stat_function(fun = dnorm, args = list(mean = ftrvar_sum$lv_med, sd = ftrvar_sum$lv_mad), col = 'red') + 
#     geom_vline(xintercept = ftrvar_sum$lv_med + 2.5 * ftrvar_sum$lv_mad, col = "blue")


# Print profiles to file --------------------------------------------------

# pdf("iprg_mq_hivar01.pdf", width = 9, height = 6)
# pdf("iprg_progenesis_hivar.pdf", width = 9, height = 6)
pdf("iprg_skyline_hivar01.pdf", width = 9, height = 6)
# pdf("bruderer_SN_hivar.pdf", width = 9, height = 6)
for (i in 1:nrow(nested_ftr)) {
    idx <- which(nested_prot2$protein == nested_ftr$protein[i])
    print(plot_profile_ftr(nested_prot2, idx, nested_ftr$data[[i]]$feature, yrange = c(10, 35)))
}
dev.off()

