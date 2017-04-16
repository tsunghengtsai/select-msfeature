# Load required packages --------------------------------------------------
library(tidyverse)
library(stringr)
library(broom)
library(survival)


# Load datasets -----------------------------------------------------------
# DDA datasets
load("data/DDA/iPRG_MaxQuant_20160919/iprg.maxquant.RData")
df_mss <- iprg.maxquant

# DIA datasets
# load("data/DIA/Bruderer_Spectronaut_20160919/bruderer.spectronaut.RData")
# df_mss <- bruderer.spectronaut


# Utility functions -------------------------------------------------------
# msstats input format to tbl to be used for feature selection
# [TODO]: truncation, censoring...
mss2ftr <- function(df_msstats) {
    df_ftr <- tbl_df(df_msstats) %>%
        select(protein = ProteinName,
               peptide = PeptideSequence,
               z_ms1 = PrecursorCharge,
               z_ms2 = ProductCharge,
               fragment = FragmentIon,
               isolab = IsotopeLabelType,
               condition = Condition,
               biorep = BioReplicate,
               run = Run,
               intensity = Intensity) %>% 
        mutate(fragment = gsub(" ", "", as.character(fragment)), 
               feature = paste(peptide, z_ms1, fragment, z_ms2, sep = "_"), 
               is_obs = !is.na(intensity), 
               log2inty = ifelse(is_obs, ifelse(intensity < 1, 0, log2(intensity)), NA)) %>% 
        mutate(protein = as.character(protein), feature = as.character(feature), run = as.character(run))
    
    return(df_ftr)
}


# convert back to the MSstats format
ftr2mss <- function(df_ftr) {
    df_msstats <- df_ftr %>% 
        mutate(protein = factor(protein), run = factor(run), fragment = factor(fragment)) %>% 
        select(ProteinName = protein, 
               PeptideSequence = peptide, 
               PrecursorCharge = z_ms1, 
               ProductCharge = z_ms2, 
               FragmentIon = fragment, 
               IsotopeLabelType = isolab, 
               Condition = condition, 
               BioReplicate = biorep, 
               Run = run, 
               Intensity = intensity, 
               IsOutlier = is_outlier) %>% 
        as.data.frame()
    
    return(df_msstats)
}

# Two-way analysis with median polish
medpolish_df <- function(df_prot) {
    # Required fields of df_prot: feature, run, log2inty
    inty_wide <- df_prot %>% select(feature, run, log2inty) %>% spread(feature, log2inty)
    inty_mat <- data.matrix(inty_wide[, -1])
    rownames(inty_mat) <- inty_wide$run
    mp_fit <- medpolish(inty_mat, na.rm = TRUE, trace.iter = FALSE)

    return(mp_fit)
}

# Generate run-level summarization from TMP fit
runsum_mp <- function(mp_fit) {
    df_sum <- tibble(
        run = names(mp_fit$row), 
        log2inty = mp_fit$overall + mp_fit$row
    )
    
    return(df_sum)
}

# Extract residuals of TMP fit
residuals_mp <- function(mp_fit) {
    df_res <- as.data.frame(mp_fit$residuals) %>% as_tibble() %>%
        rownames_to_column(var = "run") %>%
        gather(feature, log2res, -run)
    
    return(df_res)
}


# Degrees of freedom of error variance in two-way analysis
# It's only valid for linear model (as approximation for median polish)
approx_df2 <- function(df_prot) {
    df_prot <- df_prot %>% group_by(run) %>% filter(any(is_obs)) %>% ungroup()
    dfree_err <- (n_distinct(df_prot$run) - 1) * (n_distinct(df_prot$feature) - 1)
    nb_miss <- sum(!df_prot$is_obs)
    
    return(as.integer(dfree_err - nb_miss))
}


# Data wrangling ----------------------------------------------------------
# Convert to the working format
df_allftr <- mss2ftr(df_mss)

# Remove proteins with no valid measurements
df_allftr <- df_allftr %>% group_by(protein) %>% 
    filter(any(is_obs)) %>% ungroup()


# Residuals around TMP estimates ------------------------------------------
# nested_prot <- df_allftr %>% group_by(protein) %>% nest()

nested_prot <- df_allftr %>% 
    select(protein, run, feature, log2inty, is_obs) %>% 
    group_by(protein) %>% 
    nest()

nested_prot <- nested_prot %>% 
    mutate(
        mp_fit = map(data, medpolish_df), 
        mp_resid = map(mp_fit, residuals_mp), 
        mad_res = map_dbl(mp_resid, ~ mad(.$log2res, na.rm = T)), 
        nb_obs = map_int(mp_resid, ~ sum(!is.na(.$log2res))), 
        nb_ftr = map_int(data, ~ n_distinct(.$feature)), 
        dfree = map_int(data, approx_df2)
    )


# Shrinking variance estimates with limma ---------------------------------
prot_res <- nested_prot %>% 
    select(protein, mad_res, dfree) %>% 
    mutate(var_res = mad_res ^ 2)

prot_res2 <- prot_res %>% filter(mad_res != 0, dfree > 1)
em_fit <- limma::squeezeVar(var = prot_res2$var_res, df = prot_res2$dfree)
prot_res2 <- prot_res2 %>% mutate(var_eb = em_fit$var.post)


# variance estimates before/after shrinkage -------------------------------
prot_res2 %>% ggplot(aes(var_res, var_eb, colour = dfree)) + 
    geom_point() + geom_abline(slope = 1)

brk_obs <- quantile(prot_res2$dfree)
prot_res2 %>% ggplot(aes(sqrt(var_res), sqrt(var_eb), colour = dfree)) + 
    geom_point() + geom_abline(slope = 1) + 
    scale_color_gradient( trans = "log", breaks = brk_obs, labels = brk_obs)

prot_res2 %>% mutate(d_mad = sqrt(var_eb) - sqrt(var_res)) %>% 
    ggplot(aes(sqrt(var_res), d_mad, colour = dfree)) + 
    geom_point() + geom_hline(yintercept = 0) + 
    scale_color_gradient( trans = "log", breaks = brk_obs, labels = brk_obs)


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



# Different residuals -----------------------------------------------------
nested_prot <- df_allftr %>% group_by(protein) %>% nest()
nested_prot <- nested_prot %>% mutate(subdata = data)

for (i in 1:nrow(nested_prot)) {
    one_prot <- nested_prot$data[[i]]
    # Feature selection based on missing values
    obs_ftrprot <- one_prot %>% 
        select(ftr, run, log2inty) %>% 
        complete(ftr, run) %>% 
        mutate(is_ncobs = !is.na(log2inty)) %>% 
        select(-log2inty)
    mx_obs <- obs_ftrprot %>% spread(ftr, is_ncobs) %>% select(-run) %>% as.matrix()
    idx_slcftr <- select_ftrs_set(mx_obs, nb_tolmis = 0)  # index of selected features
    nb_cvr <- sum(rowSums(mx_obs) != 0)
    nb_cvr_ftr <- colSums(mx_obs)
    idx_cvrftr <- which(nb_cvr_ftr / nb_cvr >= 0.5)  # index of features with > 50% coverage
    ftr_slc <- colnames(mx_obs)[idx_slcftr]
    ftr_cvr <- colnames(mx_obs)[idx_cvrftr]
    one_prot <- one_prot %>% mutate(is_slc = ftr %in% ftr_slc, is_cvr = ftr %in% ftr_cvr)
    nested_prot$data[[i]] <- one_prot
    # Create input for AFT
    one_prot_slc <- one_prot %>% 
        filter(is_slc | is_cvr) %>% 
        group_by(ftr) %>% 
        mutate(log2inty_lim = min(log2inty[is_obs == 1], na.rm = TRUE)) %>% ungroup() %>% 
        mutate(log2inty_aft = ifelse(is_ncsr == 1, log2inty, 0.99 * log2inty_lim)) %>% 
        mutate(log2inty_na = ifelse(is_ncsr == 1, log2inty, NA)) %>% 
        select(-log2inty_lim)
    if (length(unique(one_prot_slc$ftr)) == 1) {
        one_prot_slc <- one_prot_slc %>% 
            mutate(log2inty_res = 0, res_med = 0, res_mad = 0, res_bs = 0)
    } else {
        mdl_aft <- survreg(Surv(log2inty_aft, is_ncsr, type = "left") ~ run + ftr,
                           data = one_prot_slc, dist = "gaussian")
        # pred <- predict(mdl_aft, newdata = one_prot_slc, type = "response")
        # one_prot_slc <- one_prot_slc %>%
        #     mutate(log2inty_pred = pred,
        #            log2inty_aft = ifelse(is_ncsr == 0, log2inty_pred, log2inty)) %>%
        #     select(-log2inty_pred)
        pred <- augment(mdl_aft) %>% .$.fitted
        one_prot_slc <- one_prot_slc %>%
            mutate(log2inty_aft = ifelse(is_ncsr == 0, pred, one_prot_slc$log2inty))
        
        # Residuals around lm fit
        mdl_lm <- lm(log2inty_aft ~ ftr + run, data = one_prot_slc)
        one_prot_slc <- one_prot_slc %>% mutate(log2inty_lmres = residuals(mdl_lm))
        # Spread to the wide format for TMP estimation
        wd_tmp <- one_prot_slc %>% select(ftr, run, log2inty_aft) %>% spread(ftr, log2inty_aft)
        run_name <- wd_tmp$run
        dt_tmp <- data.matrix(wd_tmp[, -1, drop = FALSE])
        dt_tmp[is.na(dt_tmp)] <- 0  # this is probably redundant
        # TMP estimate and residuals
        mp_tmp <- medpolish(dt_tmp, na.rm = TRUE, trace.iter = FALSE)  # should be no NA
        one_prot_slc <- tbl_df(mp_tmp$residuals) %>% mutate(run = run_name) %>% 
            gather(ftr, log2inty_res, -run) %>% 
            right_join(one_prot_slc, by = c("run" = "run", "ftr" = "ftr"))
        # Spread to the wide format for TMP estimation
        wd_tmp_na <- one_prot_slc %>% select(ftr, run, log2inty_na) %>% spread(ftr, log2inty_na)
        run_name_na <- wd_tmp_na$run
        dt_tmp_na <- data.matrix(wd_tmp_na[, -1, drop = FALSE])
        mp_tmp_na <- medpolish(dt_tmp_na, na.rm = TRUE, trace.iter = FALSE)
        one_prot_slc <- tbl_df(mp_tmp_na$residuals) %>% mutate(run = run_name_na) %>% 
            gather(ftr, log2inty_resna, -run) %>% 
            right_join(one_prot_slc, by = c("run" = "run", "ftr" = "ftr"))
        # Robust distance and outlier detection
        one_prot_slc <- one_prot_slc %>% 
            mutate(
                res_med = median(log2inty_res[is_obs == 1]), 
                res_mad = mad(log2inty_res[is_obs == 1]), 
                res_bs = log2inty_res / res_mad, 
                resna_med = median(log2inty_resna[is_obs == 1]), 
                resna_mad = mad(log2inty_resna[is_obs == 1])
            )
    }
    nested_prot$subdata[[i]] <- one_prot_slc
    if (i %% 100 == 0) cat(paste(i, " ")); flush.console()
}

df_allftr <- nested_prot %>% unnest(data)
df_subftr <- nested_prot %>% unnest(subdata)
# rm(nested_prot)



# Distribution of residuals -----------------------------------------------

df_subftr %>% 
    ggplot(aes(log2inty_res)) + 
    geom_histogram(binwidth = 0.01) + 
    facet_wrap(~ run) + 
    coord_cartesian(xlim = c(-1.5, 1.5))

df_subftr %>% 
    ggplot(aes(log2inty_lmres)) + 
    geom_histogram(binwidth = 0.01) + 
    facet_wrap(~ run) + 
    coord_cartesian(xlim = c(-1.5, 1.5))

df_subftr %>% distinct(protein, res_mad) %>% ggplot(aes(res_mad)) + geom_histogram(binwidth = 0.01)

df_subftr %>% filter(res_mad != 0) %>% 
    group_by(protein) %>% 
    summarise(nb_obs = sum(is_obs)) %>% 
    arrange(nb_obs) %>% 
    filter(nb_obs > 24)

df_subftr %>% filter(res_mad != 0) %>% 
    group_by(protein) %>% 
    filter(sum(is_obs) > 40) %>%
    ungroup() %>% 
    distinct(protein, res_mad) %>% 
    ggplot(aes(res_mad)) + geom_histogram(binwidth = 0.005) + coord_cartesian(xlim = c(0, 1))

df_subftr %>% filter(res_mad != 0) %>% 
    group_by(protein) %>% 
    filter(sum(is_obs) > 12) %>%
    ungroup() %>% 
    distinct(protein, res_mad) %>% 
    ggplot(aes(res_mad)) + geom_histogram(binwidth = 0.005) + coord_cartesian(xlim = c(0, 1))

df_subftr %>% distinct(protein, res_mad) %>% filter(res_mad == 0)

# Protein with only 1 feature (can't do TMP, MAD assigned as 0)
df_subftr %>% filter(res_mad == 0) %>% distinct(protein, ftr) %>% 
    count(protein) %>% filter(n != 1)

df_subftr %>% distinct(protein, res_mad) %>% filter(res_mad > 1)


nested_prot %>% filter(protein == "P53192") %>% unnest(subdata) %>% select(run, ftr, log2inty_res, log2inty_aft)
nested_prot %>% filter(protein == "Q92325") %>% unnest(subdata) %>% select(run, ftr, log2inty_res, log2inty_aft)


df_obsmad <- df_subftr %>% 
    filter(res_mad != 0) %>% 
    group_by(protein) %>% 
    summarise(res_mad = mean(res_mad), nb_obs = sum(is_obs))

ggplot(df_obsmad, aes(nb_obs, res_mad)) + 
    geom_point() + 
    geom_smooth() + 
    scale_x_log10()

df_obsmad %>% arrange(res_mad) %>% mutate(idx = seq_along(protein)) %>% 
    ggplot(aes(idx, res_mad, colour = log2(nb_obs))) + 
    geom_point() + 
    geom_line()


# Feature selection tmp residuals --------------------------------------
df_subftr <- df_subftr %>% 
    mutate(is_outlier = if_else(is_obs == 1, abs(log2inty_res) > 3 * res_mad, NA)) %>% 
    mutate(is_outlierna = if_else(is_obs == 1, abs(log2inty_resna) > 3 * resna_mad, NA))

df_subftr %>% 
    filter(is_outlier) %>% 
    select(protein, ftr, run, is_obs, log2inty_aft, log2inty_res, res_bs)

df_subftr %>% 
    filter(is_obs == 1) %>% 
    filter(res_bs != 0) %>%
    ggplot(aes(log2inty_res)) + 
    geom_histogram(binwidth = 0.01) + 
    facet_wrap(~ run) + 
    coord_cartesian(xlim = c(-1.5, 1.5))

df_subftr %>% 
    filter(is_obs == 1) %>% 
    filter(res_bs != 0) %>% 
    filter(abs(res_bs) <= 3 * bs_sd) %>% 
    ggplot(aes(res_bs)) + 
    geom_histogram(binwidth = 0.01)


# Visualization -----------------------------------------------------------
## assign run id for uncluttered visualization
uniq_run <- as.character(unique(df_subftr$run))
df_subftr$runid <- as.integer(plyr::mapvalues(df_subftr$run, from = uniq_run, to = 1:length(uniq_run)))

# With outlier
df_subftr %>% 
    filter(is_obs == 1) %>% 
    filter(res_bs != 0) %>% 
    filter(abs(res_bs) > 20 * sd(res_bs)) %>% 
    select(protein, ftr, run, is_obs, log2inty_aft, log2inty_res, res_bs) %>% 
    distinct(protein) %>% .$protein

# test --------------------------------------------------------------------
# pp <- "Q12432"
# pp <- "Q08004"
# pp <- "Q02159"
# pp <- "P53256"
# pp <- "P39729"
# pp <- "P32588"
# pp <- "P07806"
pp <- "P00635"

# pp <- "Q08004"
# pp <- "P32588"
# pp <- "P07806"
# pp <- "Q92325"
# pp <- "P53192"

dta_prot <- df_subftr %>% filter(protein == pp) %>% 
    mutate(out_cat = ifelse(!is.na(is_outlier) & is_outlier, "outlier", "non-outlier")) %>% 
    mutate(outna_cat = ifelse(!is.na(is_outlierna) & is_outlierna, "outlier", "non-outlier")) %>% 
    group_by(ftr) %>% 
    mutate(wout = any(out_cat == "outlier")) %>% 
    mutate(woutna = any(outna_cat == "outlier")) %>% 
    ungroup()
tt <- str_c(pp, round(unique(dta_prot$res_mad), 3), round(unique(dta_prot$resna_mad), 3), sep = ": ")
ggplot(dta_prot, aes(runid, log2inty_aft, group = ftr, colour = peptide, alpha = is_obs)) + 
    geom_point(aes(shape = out_cat), size = 4) +
    geom_line() + 
    scale_alpha(range = c(0.2, 0.75)) + 
    labs(x = "run", y = "log-intensity") + 
    coord_cartesian(ylim = c(10, 30)) + 
    facet_wrap(~ wout) +
    theme(legend.position = "none") + 
    ggtitle(tt)

