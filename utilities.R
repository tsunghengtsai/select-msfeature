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


# Wrapper to fit robust linear model for one protein
fit_prot_huber <- function(df_prot) {
    fit <- rlm(log2inty ~ run + feature, data = df_prot, scale.est = "Huber")
    
    return(fit)
}


# Calculate feature variance with data from broom::augment
calc_fvar <- function(augmented_data, s_resid, rm_olr = FALSE) {
    v_resid <- s_resid ^ 2
    if (rm_olr) {
        augmented_data <- augmented_data %>% 
            mutate(is_olr = abs(.resid / s_resid) > 3) %>% 
            filter(!is_olr)
    }
    varfeature <- augmented_data %>% 
        # mutate(y_samp = sample(log2inty)) %>% 
        mutate(resid_null = log2inty - mean(log2inty)) %>% 
        group_by(feature) %>%
        summarise(
            nb_run = n(), 
            svar_feature = sum(.resid ^ 2) / (nb_run - 1) / v_resid, 
            svar_ref = sum(resid_null ^ 2) / (nb_run - 1) / v_resid
            # svar_ref = var(y_samp) / v_resid
        ) %>% 
        select(feature, svar_feature, svar_ref)
    
    return(varfeature)
}


# Calculate feature variance with different ways generating ref distributions
calc_fvars <- function(augmented_data, s_resid, rm_olr = FALSE) {
    if (rm_olr) {
        augmented_data <- augmented_data %>% 
            mutate(is_olr = abs(.resid / s_resid) > 3) %>% 
            filter(!is_olr)
    }
    varfeature <- augmented_data %>% 
        mutate(y_samp = sample(log2inty)) %>% 
        mutate(resid_samp = y_samp - .fitted) %>% 
        mutate(resid_samp2 = log2inty - mean(log2inty)) %>% 
        group_by(feature) %>%
        summarise(
            nb_run = n(), 
            var_feature = sum(.resid ^ 2) / (nb_run - 1),
            var_ref = var(y_samp), 
            var_ref2 = sum(resid_samp2 ^ 2) / (nb_run - 1)
        ) %>%
        mutate(svar_feature = var_feature / s_resid ^ 2, svar_ref = var_ref / s_resid ^ 2, 
               svar_ref2 = var_ref2 / s_resid ^ 2) %>% 
        select(feature, svar_feature, svar_ref, svar_ref2)
    
    # summarise(
    #     nb_run = n(), 
    #     var_feature = sum(.resid ^ 2) / (nb_run - 1),
    #     var_ref = var(y_samp), 
    #     var_ref2 = sum(resid_samp2 ^ 2) / (nb_run - 1), 
    #     var_ref3 = var(log2inty), 
    #     var_ref4 = sum(resid_samp ^ 2) / (nb_run - 1)
    # ) %>%
    # mutate(svar_feature = var_feature / s_resid ^ 2, svar_ref = var_ref / s_resid ^ 2,
    #        svar_ref2 = var_ref2 / s_resid ^ 2, svar_ref3 = var_ref3 / s_resid ^ 2,
    #        svar_ref4 = var_ref4 / s_resid ^ 2) %>%
    # select(feature, svar_feature, svar_ref, svar_ref2, svar_ref3, svar_ref4)
    
    return(varfeature)
}



# Flag outlier with data from broom::augment
flag_outlier <- function(augmented_data, s_resid, tol = 3, keep_run = FALSE) {
    outlier <- augmented_data %>% 
        mutate(is_olr = abs(.resid / s_resid) > tol) %>% 
        select(run, feature, is_olr)
    # To keep runs from being completely removed
    if (keep_run) {
        uncovered_run <- outlier %>% 
            group_by(run) %>% 
            filter(all(is_olr)) %>% 
            ungroup() %>% 
            distinct(run)
        # Handle uncovered runs below
        if (nrow(uncovered_run) > 0) {
            outlier <- bind_rows(
                outlier %>% anti_join(uncovered_run), 
                outlier %>% semi_join(uncovered_run) %>% mutate(is_olr = FALSE)
            )
        }
    }
    
    return(outlier)
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

# Plot annotated profiles from data frame per protein
plot_profile_fltr <- function(oneprot, protein_name, show_olr = FALSE, yrange = c(10, 35), labeled = FALSE) {
    if (show_olr) {
        gprofile <- oneprot %>%
            ggplot(aes(run, log2inty, color = peptide, group = feature, shape = olr))
    } else {
        gprofile <- oneprot %>%
            ggplot(aes(run, log2inty, color = peptide, group = feature))
    }
    if (labeled) {
        gprofile <- gprofile + facet_grid(isolab ~ noise, drop = F)
    } else {
        gprofile <- gprofile + facet_grid(. ~ noise, drop = F)
    }
    
    gprofile + 
        geom_point(size = 3, alpha = 0.5) + 
        geom_line() +
        ggtitle(protein_name) +
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


# Modify bad-quality features & peaks
modify_badftrpk <- function(df_all, df_rmftr = NULL, df_rmpk = NULL, cen_val = NA, rm_olr = FALSE) {
    
    if (is.null(df_rmftr)) {
        message("df_rmftr is missing - no information about features for removal")
    }
    if (is.null(df_rmpk)) {
        message("df_rmpk is missing - no information about peaks for removal")
    }
    
    df_all <- as_tibble(df_all) %>% 
        modify_if(is.factor, as.character) %>% 
        select(protein = PROTEIN, feature = FEATURE, label = LABEL, 
               log2inty = ABUNDANCE, is_censored = censored, run = RUN, 
               group = GROUP, subject = SUBJECT, orig_run = originalRUN)
    
    # Both df_rmftr & df_rmpk are missing - return the reformated data
    if (is.null(df_rmftr) && is.null(df_rmpk)) {
        df_all$log2inty[df_all$is_censored] <- cen_val
        return(df_all)
    }
    
    # Remove noisy features
    if (!is.null(df_rmftr) && nrow(df_rmftr) > 0) {
        df_rmftr <- df_rmftr %>% 
            mutate(ftr_lab = str_c(feature_rm, label, sep = "_"))
        df_sub <- df_all %>% 
            mutate(ftr_lab = str_c(feature, label, sep = "_")) %>% 
            filter(!(ftr_lab %in% df_rmftr$ftr_lab)) %>% 
            select(-ftr_lab)
    } else {
        df_sub <- df_all
    }
    
    # Handle outliers
    if (!is.null(df_rmpk) && nrow(df_rmpk) > 0) {
        df_rmpk <- df_rmpk %>% 
            mutate(peak_lab = str_c(run, feature, label, sep = "_"))
        df_sub <- df_sub %>% 
            mutate(peak_lab = str_c(orig_run, feature, label, sep = "_"))
        
        if (rm_olr) {
            df_sub <- df_sub %>% 
                filter(!(peak_lab %in% df_rmpk$peak_lab)) %>% 
                select(-peak_lab)
        } else {
            # Censored the outliers of the endogenous (may need revisions...)
            df_sub$is_censored[df_sub$label == "L" & df_sub$peak_lab %in% df_rmpk$peak_lab] <- TRUE
            # Replace outliers of the labeled with NA
            df_sub$log2inty[df_sub$label == "H" & df_sub$peak_lab %in% df_rmpk$peak_lab] <- NA
            
            df_sub <- df_sub %>% select(-peak_lab)
        }
    }
    df_sub$log2inty[df_sub$is_censored] <- cen_val
    
    return(df_sub)
}


# Generate subplot summarization & intermediate results (following modify_badftrpk)
gen_subplot <- function(df_mss, censoredInt = NULL, MBimpute = TRUE, message = TRUE) {
    
    if (is.null(censoredInt)) {
        message("censoredInt is not provided - set as 'NA'")
        censoredInt <- "NA"
    } else if (length(censoredInt) > 1) {
        error("censoredInt should be either '0' or 'NA'")
    } else if (!(censoredInt %in% c("NA", "0"))) {
        error("Set censoredInt as '0' or 'NA'")
    }
    
    # Format
    df_mss <- as_tibble(df_mss) %>% 
        modify_if(is.factor, as.character) %>% 
        mutate(
            ftr_lab = paste(feature, label, sep = "_"), 
            run_lab = paste(run, label, sep = "_")
        )
    
    # Is the dataset label-based (or label-free)?
    is_labeled <- n_distinct(df_mss$label) == 2
    
    # Replace censored values with censoredInt
    if (MBimpute) {
        if (censoredInt == "0") {
            df_mss$log2inty[df_mss$is_censored] <- 0
        } else if (censoredInt == "NA") {
            df_mss$log2inty[df_mss$is_censored] <- NA
        }
        df_mss <- df_mss %>% 
            mutate(ind_obs = ifelse(is_censored, 0, 1))
    }
    df_mss <- df_mss %>% 
        mutate(is_missed = is.na(log2inty))
    
    # Warning when a measurement is <= 0 but uncensored
    if (df_mss %>% filter(log2inty <= 0, !is_censored) %>% nrow() > 0) {
        message("There are abundance values less than or equal to 0 but uncensored!")
    }
    
    # Create a nested data frame wrt protein
    nested <- df_mss %>% nest(-protein)
    nb_prot <- nrow(nested)
    
    # Empty lists to save intermediate results (assuming MBimpute)
    if (is_labeled) {
        res_aft <- res_subplot <- res_labeled <- res_adjlabeled <- vector("list", nb_prot)
    } else {
        res_aft <- res_subplot <- vector("list", nb_prot)
    }
    
    for (i in 1:nb_prot) {
        
        str_prot <- paste(nested$protein[i], "(", i, "of", nb_prot, ")")
        if (message) {
            message(paste("Getting the summarization by Tukey's median polish for protein", str_prot))
        }
        
        oneprot <- nested$data[[i]]
        # Check if the protein has valid peaks for summarization
        valid_pks <- oneprot %>% 
            filter(label == "L", !is_missed, !is_censored)
        if (nrow(valid_pks) == 0) {
            message(paste("Can't summarize for", str_prot, "because all measurements are NAs"))
            next()
        }
        
        # Remove features with only one measurement
        cnt_feature <- valid_pks %>% 
            count(feature)
        if (all(cnt_feature$n == 1)) {
            message(paste("Can't summarize for", str_prot, "because features have only one measurement across MS runs"))
            next()
        }
        # [TODO]: this would also remove labeled peptides - to be confirmed
        valid_ftrs <- cnt_feature %>% 
            filter(n > 1) %>% 
            .$feature
        oneprot <- oneprot %>% filter(feature %in% valid_ftrs)
        
        # Remove uncovered runs (all endogenous are either missing or censored)
        covered <- valid_pks %>% 
            filter(feature %in% valid_ftrs) %>% 
            distinct(run) %>% 
            .$run
        oneprot <- oneprot %>% filter(run %in% covered)
        
        # Impute censored values
        if (MBimpute) {
            if (is_labeled) {
                oneprot_h <- oneprot %>% filter(label == "H")
                oneprot_l <- oneprot %>% filter(label == "L")
            } else {
                oneprot_l <- oneprot
            }
            if (any(oneprot_l$is_censored)) {
                # Replace censored values by feature-specific thresholds
                oneprot_l <- oneprot_l %>% 
                    group_by(ftr_lab) %>% 
                    mutate(cen_thresh = 0.99 * min(log2inty[!(is_censored | is_missed)])) %>% 
                    ungroup() %>% 
                    mutate(log2inty = ifelse(is_censored, cen_thresh, log2inty)) %>% 
                    select(-cen_thresh)
                
                # AFT model
                valid_aft <- oneprot_l %>% filter(!is_missed)
                insufDF <- nrow(valid_aft) < (n_distinct(valid_aft$feature) + n_distinct(valid_aft$run) - 1)
                
                set.seed(100)
                if (n_distinct(valid_aft$feature) == 1) {
                    aftfit <- survreg(Surv(log2inty, ind_obs, type = "left") ~ run, data = oneprot_l, dist = "gaussian")
                } else {
                    if (insufDF) {
                        aftfit <- survreg(Surv(log2inty, ind_obs, type = "left") ~ run, data = oneprot_l, dist = "gaussian")
                    } else {
                        aftfit <- survreg(Surv(log2inty, ind_obs, type = "left") ~ feature + run, data = oneprot_l, dist = "gaussian")
                    }
                }
                # Replace censored values with predicted values by AFT
                oneprot_l <- oneprot_l %>% 
                    mutate(pred = predict(aftfit, newdata = oneprot_l, type = "response")) %>% 
                    mutate(log2inty = ifelse(is_censored, pred, log2inty))
                
                res_aft[[i]] <- oneprot_l %>%
                    mutate(protein = nested$protein[i])
                
                if (is_labeled) {
                    # Merge with labeled peptide
                    oneprot_h$pred <- NA
                    oneprot <- rbind(oneprot_h, oneprot_l)
                } else {
                    oneprot <- oneprot_l
                }
                
            } else {
                # No censored values
                res_aft[[i]] <- oneprot_l %>%
                    mutate(pred = NA, protein = nested$protein[i])
                oneprot$pred <- NA
            }
        }
        # Remove NA in abundance prior to TMP (likely introduced by AFT?)
        oneprot <- oneprot %>% filter(!is.na(log2inty))
        
        # Tukey's median polish
        if (n_distinct(oneprot$feature) > 1) { 
            # More than 1 features are required for TMP
            if (!is_labeled) {
                inty_wide <- oneprot %>% 
                    select(feature, run, log2inty) %>% 
                    spread(feature, log2inty)
                inty_mat <- data.matrix(inty_wide[, -1])
                mp_out <- medpolish(inty_mat, na.rm = TRUE, trace.iter = FALSE)
                
                # Subplot result
                res_subplot[[i]] <- tibble(
                    protein = nested$protein[i], 
                    log2inty = mp_out$overall + mp_out$row, 
                    run = inty_wide$run
                )
            } else { 
                # Labeled - additional adjustment based on the reference
                inty_wide <- oneprot %>% 
                    select(feature, run_lab, log2inty) %>% 
                    spread(feature, log2inty)
                inty_mat <- data.matrix(inty_wide[, -1])
                mp_out <- medpolish(inty_mat, na.rm = TRUE, trace.iter = FALSE)
                
                df_labeled <- tibble(
                    protein = nested$protein[i], 
                    log2inty = mp_out$overall + mp_out$row, 
                    run_lab = inty_wide$run_lab
                ) %>% 
                    mutate(
                        run = str_sub(run_lab, 1, -3), 
                        label = str_sub(run_lab, -1, -1)
                    ) %>% 
                    select(-run_lab)
                
                # Adjustment wrt the labeled peptides
                adj_h <- df_labeled %>% 
                    filter(label == "H") %>% 
                    mutate(log2inty_adj = median(log2inty, na.rm = TRUE) - log2inty) %>% 
                    select(run, log2inty_adj)
                df_adjlabeled <- df_labeled %>% 
                    left_join(adj_h) %>% 
                    mutate(log2inty = ifelse(!is.na(log2inty_adj), log2inty + log2inty_adj, log2inty)) %>% 
                    select(-log2inty_adj)
                
                res_labeled[[i]] <- df_labeled
                res_adjlabeled[[i]] <- df_adjlabeled
                
                # Adjusted subplot results for the endogenous
                res_subplot[[i]] <- df_adjlabeled %>% 
                    filter(label == "L") %>% 
                    select(-label)
            }
        } else { 
            # single feature 
            if (!is_labeled) {
                res_subplot[[i]] <- oneprot %>% 
                    filter(!is.na(log2inty)) %>% 
                    mutate(protein = nested$protein[i]) %>% 
                    select(protein, log2inty, run)
            } else {
                # Labeled - additional adjustment
                df_labeled <- oneprot %>% 
                    filter(!is.na(log2inty)) %>% 
                    mutate(protein = nested$protein[i]) %>% 
                    select(protein, log2inty, run, label)
                adj_h <- df_labeled %>% 
                    filter(label == "H") %>% 
                    mutate(log2inty_adj = median(log2inty, na.rm = TRUE) - log2inty) %>% 
                    select(run, log2inty_adj)
                df_adjlabeled <- df_labeled %>% 
                    left_join(adj_h) %>% 
                    mutate(log2inty = ifelse(!is.na(log2inty_adj), log2inty + log2inty_adj, log2inty)) %>% 
                    select(-log2inty_adj)
                
                res_labeled[[i]] <- df_labeled
                res_adjlabeled[[i]] <- df_adjlabeled
                
                # Adjusted subplot results for the endogenous
                res_subplot[[i]] <- df_adjlabeled %>% 
                    filter(label == "L") %>% 
                    select(-label)
            }
        }
    }
    
    # Output
    if (is_labeled) {
        output <- list(
            res_aft = bind_rows(res_aft), 
            res_subplot = bind_rows(res_subplot), 
            res_labeled = bind_rows(res_labeled), 
            res_adjlabeled = bind_rows(res_adjlabeled)
        )
    } else {
        output <- list(
            res_aft = bind_rows(res_aft), 
            res_subplot = bind_rows(res_subplot)
        )
    }
    
    return(output)
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


# Calculate feature variance with data from broom::augment
calc_varfeature <- function(augmented_data, s_resid) {
    varfeature <- augmented_data %>%
        mutate(y_samp = sample(log2inty)) %>%
        mutate(resid_samp = y_samp - .fitted) %>%
        group_by(feature) %>%
        summarise(
            nb_run = n(),
            var_feature = sum(.resid ^ 2) / (nb_run - 1),
            var_ref = sum(resid_samp ^ 2) / (nb_run - 1)
        ) %>%
        mutate(svar_feature = var_feature / s_resid ^ 2, svar_ref = var_ref / s_resid ^ 2) %>%
        select(feature, svar_feature, svar_ref)
    
    return(varfeature)
}
