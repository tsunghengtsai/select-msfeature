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
        rlm_fit = map(subdata, possibly(~ rlm(log2inty ~ 0 + run + feature, data = ., scale.est = "Huber"), NULL)), 
        dfree = map_dbl(subdata, approx_df2)
    )


# L1-penalty least absolute deviation fit ---------------------------------

library(flare)

min_l1lad <- function(rlm_fit, subdata) {
    resp <- rlm_fit %>% tidy() %>% filter(str_detect(term, "run")) %>% .$estimate
    fulldata <- subdata %>% select(run, feature, log2inty) %>% complete(run, feature)
    fitted <- rlm_fit %>% augment(., newdata = fulldata) %>% 
        mutate(log2inty_fit = ifelse(is.na(log2inty), .fitted, log2inty))
    predvar <- fitted %>% select(run, feature, log2inty_fit) %>% 
        spread(feature, log2inty_fit) %>% select(-run) %>% as.matrix()
    reg <- 0.5 * sqrt(log(ncol(predvar)) / nrow(predvar))
    l1lad <- slim(X = predvar, Y = resp, method = "lq", q = 1, lambda = reg, verbose = FALSE)
    
    return(l1lad)
}

nested_prot <- nested_prot %>% 
    filter(dfree > 0) %>% 
    mutate(l1lad_fit = map2(rlm_fit, subdata, min_l1lad))

nested_prot <- nested_prot %>% 
    mutate(feature = map(l1lad_fit, ~ colnames(.$X)[.$beta == 0]))

nested_ftr <- nested_prot %>% unnest(feature) %>% nest(-protein)


# Print profiles to file --------------------------------------------------

pdf("iprg_mq_l1lad01.pdf", width = 9, height = 6)
# pdf("iprg_progenesis_l1lad01.pdf", width = 9, height = 6)
# pdf("iprg_skyline_l1lad01.pdf", width = 9, height = 6)
# pdf("bruderer_SN_l1lad01.pdf", width = 9, height = 6)
for (i in 1:nrow(nested_ftr)) {
    idx <- which(nested_prot$protein == nested_ftr$protein[i])
    print(plot_profile_ftr(nested_prot, idx, nested_ftr$data[[i]]$feature, yrange = c(10, 35)))
}
dev.off()

