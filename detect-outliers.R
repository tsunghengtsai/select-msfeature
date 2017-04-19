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


# Flag features with low and non-essential coverage
flag_lowcover <- function(df_prot) {
    mat_obsftr <- df_prot %>% select(run, feature, is_obs) %>% 
        spread(feature, is_obs) %>% 
        select(-run) %>% as.matrix()
    is_hicover <- (colSums(mat_obsftr) / nrow(mat_obsftr)) > 0.5
    if (all(is_hicover)) {
        # All features have high coverage
        ftr_hicover <- colnames(mat_obsftr)[is_hicover]
    } else {
        # Otherwise consider set cover 
        ftr_hicover <- colnames(mat_obsftr)[union(idx_coverset(mat_obsftr), which(is_hicover))]
    }
    
    return(df_prot %>% mutate(is_lowcvr = !(feature %in% ftr_hicover)))
}


# Index of columns (features) covering complete set of rows (runs)
# Matrix of logical variables as input
idx_coverset <- function(mat_obs) {
    if (ncol(mat_obs) == 1) return(1)  # only one feature
    cost <- colSums(!mat_obs)  # number of missing values as cost
    # return all fully covering features
    if (any(cost == 0)) return(which(cost == 0))
    is_complete <- FALSE
    set_ftr <- NULL
    # find features with the smallest cost (# of missing / new additions)
    # until no more additions can be found (all runs are covered)
    # NB: runs covered by selected features will be removed from mat_obs
    while (!is_complete) {
        if (nrow(mat_obs) == 1) {
            add <- c(mat_obs)
        } else {
            add <- colSums(mat_obs)
        }
        if (all(add == 0)) {
            is_complete <- TRUE
        } else {
            cost_per_add <- cost / add
            min_cost <- min(cost_per_add)
            new_ftr <- which(cost_per_add == min_cost)
            set_ftr <- c(set_ftr, new_ftr)
            if (length(new_ftr) == 1) {
                idx_add <- which(mat_obs[, new_ftr])
            } else {
                idx_add <- which(rowSums(mat_obs[, new_ftr, drop = FALSE]) != 0)
            }
            mat_obs[, new_ftr] <- FALSE  # exclude selected features from further consideration
            mat_obs <- mat_obs[-idx_add, , drop = FALSE]  # remove covered runs
            if (nrow(mat_obs) == 0) is_complete <- TRUE
        }
    }
    
    return(set_ftr)
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


# Degrees of freedom for MAD estimate in two-way analysis
# MAD is 37% as efficient as sample SD
approx_df2_mad <- function(df_prot) {
    df_prot <- df_prot %>% group_by(run) %>% filter(any(is_obs)) %>% ungroup()
    dfree_err <- (n_distinct(df_prot$run) - 1) * (n_distinct(df_prot$feature) - 1)
    nb_miss <- sum(!df_prot$is_obs)
    df2 <- dfree_err - nb_miss
    # df2_mad <- floor(max(1, df2 * 0.37))
    df2_mad <- max(1, df2 * 0.37)
    
    return(df2_mad)
}


# Select features based on the flags of outliers
# Features are flagged as high-variance if their non-outlier measurements 
# do not provide reasonable coverage (50%)
# flag_hivar <- function(df_res) {
#     # The matrix may contain NA (can't determine outlyingness for unobserved data)
#     mat_olrftr <- df_res %>% select(run, feature, is_olr) %>% 
#         spread(feature, is_olr) %>% 
#         select(-run) %>% as.matrix()
#     # Test for high-coverage features (non-outlying observed measurements account for >50%)
#     is_hicover <- (colSums(!mat_olrftr, na.rm = TRUE) / nrow(mat_olrftr)) > 0.5
#     
#     if (any(is_hicover)) {
#         mat_hicover <- mat_olrftr[, is_hicover, drop = FALSE]
#         
#         # Test if all runs covered by >1 high-coverage features
#         run_cover <- rowSums(!is.na(mat_hicover)) > 0
#         run_incover <- rowSums(!mat_hicover, na.rm = TRUE) > 0
#         if (all(run_cover) && all(run_incover)) {
#             # All covered by high-coverage features
#             set_ftr <- colnames(mat_hicover)
#         } else if (all(run_cover) && any(run_incover)) {
#             # Some runs become uncovered when removing outliers -> swap good candidates
#             set_ftr <- colnames(mat_hicover)
#             run_per_ftr <- colSums(!mat_hicover, na.rm = TRUE)
#             for (i in which(!run_incover)) {
#                 idx_cand <- which(!is.na(mat_hicover[i, ]))
#                 run_per_cand <- run_per_ftr[idx_cand]
#                 idx_maxcand <- which(run_per_cand == max(run_per_cand))
#                 mat_hicover[i, idx_cand[idx_maxcand]] <- FALSE
#             }
#         }
# 
#     }
# 
#     return(df_res %>% mutate(is_hivar = !(feature %in% set_ftr)))
# }



flag_hivar <- function(df_res) {
    # The matrix may contain NA (can't determine outlyingness for unobserved data)
    mat_olrftr <- df_res %>% select(run, feature, is_olr) %>% 
        spread(feature, is_olr) %>% 
        select(-run) %>% as.matrix()
    # Test for no-outlier ("inlier") features
    is_allin <- colSums(mat_olrftr, na.rm = TRUE) == 0
    if (any(is_allin)) {
        mat_inftr <- mat_olrftr[, is_allin, drop = FALSE]
        mat_outftr <- mat_olrftr[, !is_allin, drop = FALSE]
        # Test if all runs covered by at least one no-outlier features
        run_cover <- rowSums(!is.na(mat_inftr)) > 0
        # is_incover <- all(rowSums(!is.na(mat_inftr)) > 0)
        is_hicover <- (colSums(!mat_outftr, na.rm = TRUE) / nrow(mat_outftr)) > 0.5
        if (all(run_cover)) {
            # Test if non-outlier (& observed) measurements account for >50%
            # is_hicover <- (colSums(!mat_outftr, na.rm = TRUE) / nrow(mat_outftr)) > 0.5
            set_ftr <- c(colnames(mat_inftr), colnames(mat_outftr)[is_hicover])
        } else {
            # Acquire key players
            # idx_misrun <- which(!run_cover)
            set_ftr <- c(colnames(mat_inftr), colnames(mat_outftr)[is_hicover])
        }
    } else {
        # All features contain at least one outlier
        # (19, 20, 111, 130, 1846, 2521)
        mat_outftr <- mat_olrftr[, !is_allin, drop = FALSE]
        is_hicover <- (colSums(!mat_outftr, na.rm = TRUE) / nrow(mat_outftr)) > 0.5
        ftr_cover <- colnames(mat_outftr)[is_hicover]
    }
    
    return(df_res %>% mutate(is_hivar = !(feature %in% set_ftr)))
}



# Data wrangling ----------------------------------------------------------
# Convert to the working format
df_allftr <- mss2ftr(df_mss)

# For each protein, remove uncovered runs and undetectable features
# Subsequent processing assumes this is done
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
        nb_obs = map_int(mp_resid, ~ sum(!is.na(.$log2res))), 
        dfree = map_dbl(subdata, approx_df2_mad)
    )


# Shrinkage variance estimation with limma --------------------------------
prot_res <- nested_prot %>% 
    select(protein, mad_res, dfree) %>% 
    mutate(var_res = mad_res ^ 2)

# prot_res2 <- prot_res %>% filter(mad_res != 0, dfree > 1)
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

nested_prot <- nested_prot %>% 
    mutate(mp_resid = map(mp_resid, flag_hivar))

nested_prot %>% mutate(rmv = map_lgl(mp_resid, ~ all(.$is_hivar))) %>% filter(rmv)

# Plot profiles with the flag ---------------------------------------------
plot_annoprofile <- function(nested, idx) {
    oneprot <- nested$data[[idx]] %>% left_join(nested$mp_resid[[idx]]) %>% 
        mutate(olr = if_else(is_olr, "Yes", "No/Unsure", "No/Unsure")) %>% 
        group_by(feature) %>% 
        mutate(
            ftr_olr = if_else(is_lowcvr, "Low coverage", if_else(any(is_olr, na.rm = T), "W/ outlier", "W/O outlier")), 
            ftr_olr = factor(ftr_olr, levels = c("Low coverage", "W/ outlier", "W/O outlier"))
        ) %>% 
        ungroup()
    oneprot %>% 
        ggplot(aes(run, log2inty, color = peptide, group = feature, shape = olr, alpha = olr)) + 
        geom_point(size = 3) + 
        geom_line() + 
        scale_alpha_discrete(range = c(0.5, 1)) + 
        ggtitle(nested$protein[idx]) + 
        facet_wrap(~ ftr_olr, drop = F) + 
        coord_cartesian(ylim = c(17, 35)) + 
        theme(legend.position = "none", 
              # axis.text.x = element_text(angle = 90, hjust = 1), 
              axis.text.x = element_blank())
}


plot_annoprofile2 <- function(nested, idx) {
    oneprot <- nested$data[[idx]] %>% left_join(nested$mp_resid[[idx]]) %>% 
        mutate(olr = if_else(is_olr, "Yes", "No/Unsure", "No/Unsure")) %>% 
        mutate(ftr_olr = if_else(is_lowcvr | is_hivar, "Remove", "Select"), 
               ftr_olr = factor(ftr_olr, levels = c("Remove", "Select")))
    oneprot %>% 
        ggplot(aes(run, log2inty, color = peptide, group = feature, shape = olr, alpha = olr)) + 
        geom_point(size = 3) + 
        geom_line() + 
        scale_alpha_discrete(range = c(0.5, 1)) + 
        ggtitle(nested$protein[idx]) + 
        facet_wrap(~ ftr_olr, drop = F) + 
        coord_cartesian(ylim = c(17, 35)) + 
        theme(legend.position = "none", 
              axis.text.x = element_blank())
}

# Print profiles to file --------------------------------------------------

# pdf("profile_out_10.pdf", width = 9, height = 6)
# for (i in 901:1000) {
#     print(plot_annoprofile(nested_prot, i))
# }
# dev.off()

pdf("profile_slc_02.pdf", width = 9, height = 6)
for (i in 101:200) {
    print(plot_annoprofile2(nested_prot, i))
}
dev.off()


# Spike-in ----------------------------------------------------------------
prot_spike <- c("P44015", "P55752", "P44374", "P44983", "P44683", "P55249")

pdf("profile_out_spike.pdf", width = 9, height = 6)
for (i in which(nested_prot$protein %in% prot_spike)) {
    print(plot_annoprofile(nested_prot, i))
}
dev.off()

pdf("profile_slc_spike.pdf", width = 9, height = 6)
for (i in which(nested_prot$protein %in% prot_spike)) {
    print(plot_annoprofile2(nested_prot, i))
}
dev.off()


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

