# Processed data to tbl for feature selection
pdata2ftr <- function(pdata) {
    df_ftr <- tbl_df(pdata) %>%
        select(
            protein = PROTEIN,
            peptide = PEPTIDE,
            feature = FEATURE, 
            run = originalRUN, 
            isolab = LABEL,
            log2inty = ABUNDANCE, 
            is_censored = censored
        ) %>% 
        mutate(
            is_obs = !(is.na(log2inty) | is_censored), 
            log2inty = ifelse(is_obs, log2inty, NA), 
            protein = as.character(protein), 
            peptide = as.character(peptide), 
            feature = as.character(feature), 
            run = as.character(run)
        )
    
    return(df_ftr)
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
flag_hivar <- function(df_res) {
    # The matrix may contain NA (can't determine outlyingness for unobserved data)
    # or bc var_eb, var_reb are NA (no way to decide...)
    if (all(is.na(df_res$is_olr))) return(df_res %>% mutate(is_hivar = FALSE))
    mat_olrftr <- df_res %>% select(run, feature, is_olr) %>%
        spread(feature, is_olr) %>%
        select(-run) %>% as.matrix()
    # Test for high-coverage features (non-outlying observed measurements account for >50%)
    is_hicover <- (colSums(!mat_olrftr, na.rm = TRUE) / nrow(mat_olrftr)) > 0.5
    set_ftr <- colnames(mat_olrftr)[is_hicover]

    return(df_res %>% mutate(is_hivar = !(feature %in% set_ftr)))
}


# List features to be removed
list_rmfeature <- function(df_prot, df_res) {
    ftr_lowcvr <- df_prot %>% filter(is_lowcvr) %>% distinct(feature) %>% .$feature
    ftr_hivar <- df_res %>% filter(is_hivar) %>% distinct(feature) %>% .$feature
    
    return(c(ftr_lowcvr, ftr_hivar))
}


# Plot annotated profiles from nested data frame
plot_profile_nest <- function(nested, idx, yrange = c(10, 35)) {
    oneprot <- nested$data[[idx]] %>% left_join(nested$subdata[[idx]]) %>% 
        mutate(olr = if_else(is_olr, "Yes", "No/Unsure", "No/Unsure")) %>% 
        mutate(ftr_olr = if_else(is_lowcvr | is_hivar, "Remove", "Select"), 
               ftr_olr = factor(ftr_olr, levels = c("Remove", "Select")))
    oneprot %>% 
        ggplot(aes(run, log2inty, color = peptide, group = feature, shape = olr, alpha = olr)) + 
        geom_point(size = 3) + geom_line() + 
        scale_alpha_discrete(range = c(0.5, 1)) + 
        ggtitle(nested$protein[idx]) + 
        facet_wrap(~ ftr_olr, drop = F) + 
        coord_cartesian(ylim = yrange) + 
        theme(legend.position = "none", axis.text.x = element_blank())
}


# Plot annotated profiles from nested data frame
plot_profile_nest2 <- function(nested, idx, yrange = c(10, 35)) {
    oneprot <- nested$data[[idx]] %>% left_join(nested$subdata[[idx]]) %>% 
        mutate(olr = if_else(is_olr, "Yes", "No/Unsure", "No/Unsure")) %>% 
        group_by(feature) %>% 
        mutate(ftr_olr = if_else(all(is_lowcvr) | all(is_hivar), "Remove", 
                                 if_else(any(olr == "Yes"), "w/ outliers", "w/o outliers"))) %>% 
        ungroup() %>% 
        mutate(ftr_olr = factor(ftr_olr, levels = c("Remove", "w/ outliers", "w/o outliers")))
    oneprot %>% 
        ggplot(aes(run, log2inty, color = peptide, group = feature, shape = olr, alpha = olr)) + 
        geom_point(size = 3) + geom_line() + 
        scale_alpha_discrete(range = c(0.5, 1)) + 
        ggtitle(nested$protein[idx]) + 
        facet_wrap(~ ftr_olr, drop = F) + 
        coord_cartesian(ylim = yrange) + 
        theme(legend.position = "none", axis.text.x = element_blank())
}


# Experimental ------------------------------------------------------------

# Extract scaled residuals
extract_sres <- function(fit) {
    shat <- summary(fit)$sigma
    fit %>% 
        augment() %>% 
        mutate(scl_r = .resid / shat, std_r = scl_r / sqrt(1 - .hat)) %>% 
        select(run, feature, scl_r, std_r)
}

# Posterior probability of two-component mixture model
calc_mixprob <- function(x, mix_mu, mix_sigma, mix_lambda) {
    logodd <- dnorm(x, mix_mu[2], mix_sigma[2], log = TRUE) - 
        dnorm(x, mix_mu[1], mix_sigma[1], log = TRUE) + 
        log(mix_lambda[2]) - log(mix_lambda[1])
    
    return(1 / (1 + exp(logodd)))
}

# Scaled normal density
snorm <- function(x, mu, sigma, lambda) {
    dnorm(x, mu, sigma) * lambda
}


# Extract residuals & scaled residuals
extract_res <- function(fit) {
    shat <- summary(fit)$sigma
    fit %>% 
        augment() %>% 
        mutate(r = .resid, scl_r = .resid / shat) %>% 
        select(run, feature, r, scl_r)
}


# Being deprecated --------------------------------------------------------

# Intensity column is not normalized across runs
pdata2ftr_old <- function(pdata) {
    df_ftr <- tbl_df(pdata) %>%
        select(
            protein = PROTEIN,
            peptide = PEPTIDE,
            feature = FEATURE, 
            run = originalRUN, 
            isolab = LABEL,
            intensity = INTENSITY, 
            is_censored = censored
        ) %>% 
        mutate(
            is_obs = !(is.na(intensity) | is_censored), 
            log2inty = ifelse(is_obs, ifelse(intensity < 1, 0, log2(intensity)), NA), 
            protein = as.character(protein), 
            peptide = as.character(peptide), 
            feature = as.character(feature), 
            run = as.character(run)
        )
    
    return(df_ftr)
}


# msstats input format to tbl to be used for feature selection
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


