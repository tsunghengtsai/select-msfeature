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


# Fit with linear model for group comparison
fit_grpcomp <- function(one_protein) {
    fit <- lm(log2inty ~ 0 + group, data = one_protein)
    
    return(fit)
}


# Summary of estimated parameters
tidy_grpcomp <- function(fit) {
    params <- tidy(fit) %>% 
        mutate(group = str_replace(term, "group", "")) %>% 
        select(group, estimate, std.error)
    
    return(params)
}


# Testing for differential abundance
test_grpcomp <- function(all_protein, grp_ctrl, grp_case, grp_code = "") {
    filled <- all_protein %>% 
        filter(group %in% c(grp_ctrl, grp_case)) %>% 
        filter(!is.na(df_res)) %>% 
        mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>% 
        complete(protein, key_grp)
    # Model-based inference of log-difference
    full_ctrx <- filled %>% 
        group_by(protein) %>% 
        filter(!any(is.na(df_res))) %>% 
        group_by(protein, df_res) %>% 
        summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>% 
        ungroup()
    res <- full_ctrx %>% 
        mutate(logFC = log2fc, SE = sqrt(se2), DF = df_res, t_stat = logFC / SE, 
               p_val = 2 * pt(abs(t_stat), df = DF, lower.tail = FALSE)) %>% 
        select(protein, logFC, SE, DF, t_stat, p_val) %>% 
        mutate(ctrx = paste0(grp_code, grp_case, "-", grp_code, grp_ctrl))
    # Missing in one group
    part_ctrx <- filled %>% 
        anti_join(full_ctrx %>% select(protein)) %>% 
        filter(!is.na(df_res)) %>% 
        distinct(protein, key_grp)
    if (nrow(part_ctrx) > 0) {
        untest <- part_ctrx %>% 
            mutate(
                logFC = ifelse(key_grp == "G1", Inf, -Inf), 
                SE = NA, DF = NA, t_stat = NA, p_val = NA, 
                ctrx = paste0(grp_code, grp_case, "-", grp_code, grp_ctrl)
            ) %>% select(-key_grp)
        res <- res %>% bind_rows(untest)
    }
    
    return(res)
}


# Evaluation
sensPrecCurves <- function(test_results) {
    meths <- unique(test_results$method)
    test_meth = vector("list", length = length(meths))
    for (i in seq_along(meths)) {
        test_sub <- test_results %>% 
            filter(method == meths[i], !is.na(pvalue)) %>% 
            arrange(pvalue)
        ordered_de <- test_sub$de == "With change"
        test_meth[[i]] <- tibble(
            sens = cumsum(ordered_de) / sum(ordered_de), 
            fdr = cumsum(1 - ordered_de) / seq_along(ordered_de),
            fpr = cumsum(!ordered_de) / sum(!ordered_de)
        ) %>% mutate(method = meths[i])
    }
    
    return(bind_rows(test_meth))
}
