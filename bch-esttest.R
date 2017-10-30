
# Whole-plot modeling & testing -------------------------------------------

library(tidyverse)
library(stringr)
library(broom)

source("utilities-downstream.R")

load(paste0("fltr/", dset_name, ".rda"))
eval(parse(text = paste0("df_mss <- ", dset_name)))
df_orig <- modify_badftrpk(df_mss)

res_orig <- readRDS(paste0("subplot/", str_replace_all(dset_name, "\\.", "_"), ".rds"))
sub_orig <- res_orig$res_subplot

# Experimental design
design <- df_orig %>% distinct(group, subject, run)

# Nested data frames
nested_orig <- sub_orig %>% 
    left_join(design) %>% 
    nest(-protein)

# Modeling & parameter estimation
nested_orig <- nested_orig %>% 
    mutate(lm_fit = map(data, possibly(fit_grpcomp, NULL))) %>% 
    mutate(fail = map_lgl(lm_fit, is.null)) %>% 
    filter(!fail) %>% 
    select(-fail) %>% 
    mutate(param = map(lm_fit, tidy_grpcomp)) %>% 
    mutate(df_res = map_dbl(lm_fit, df.residual))

param_orig <- nested_orig %>% 
    unnest(param)

vers <- c("v2", "v3", "v4", "v5", "v6", "v7")
param_fltrs <- vector("list", length(vers))
for (i in seq_along(vers)) {
    meth_ver <- vers[i]
    res_fltr <- readRDS(paste0("subplot/", str_replace_all(dset_name, "\\.", "_"), "_m5", meth_ver, ".rds"))
    sub_fltr <- res_fltr$res_subplot
    
    nested_fltr <- sub_fltr %>% 
        left_join(design) %>% 
        nest(-protein)
    
    nested_fltr <- nested_fltr %>% 
        mutate(lm_fit = map(data, possibly(fit_grpcomp, NULL))) %>% 
        mutate(fail = map_lgl(lm_fit, is.null)) %>% 
        filter(!fail) %>% 
        select(-fail) %>% 
        mutate(param = map(lm_fit, tidy_grpcomp)) %>% 
        mutate(df_res = map_dbl(lm_fit, df.residual))
    
    param_fltrs[[i]] <- nested_fltr %>% 
        unnest(param)
}


# Model-based testing
if (str_detect(dset_name, "choi")) {
    cases = c("1", "1", "1", "1", "2", "2", "2", "3", "3", "4")
    controls = c("2", "3", "4", "5", "3", "4", "5", "4", "5", "5")
    cond_code <- "M"
} else if (str_detect(dset_name, "iprg")) {
    controls = c("1", "1", "1", "2", "2", "3")
    cases = c("2", "3", "4", "3", "4", "4")
    cond_code <- "C"
} else if (str_detect(dset_name, "navarro")) {
    controls = c("2")
    cases = c("1")
    cond_code <- ""
}

test_orig <- vector("list", length = length(cases))
test_fltr <- vector("list", length = length(cases))
for (i in seq_along(cases)) {
    grp_ctrl <- controls[i]
    grp_case <- cases[i]
    test_orig[[i]] <- test_grpcomp(param_orig, grp_ctrl, grp_case, cond_code) %>% 
        mutate(method = "All")
    
    test_vers <- vector("list", length = length(vers))
    for (j in seq_along(vers)) {
        test_vers[[j]] <- test_grpcomp(param_fltrs[[j]], grp_ctrl, grp_case, cond_code) %>% 
            mutate(method = paste0("m5", vers[j]))
    }
    test_fltr[[i]] <- bind_rows(test_vers)
}

test_res <- bind_rows(
    bind_rows(test_orig), 
    bind_rows(test_fltr)
)

test_res <- test_res %>% 
    group_by(method) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
    mutate(p_adj = ifelse(logFC %in% c(Inf, -Inf), 0, p_adj)) %>% 
    ungroup() %>% 
    mutate(pvalue = p_val, adj.pvalue = p_adj)

saveRDS(test_res, paste0("esttest/", str_replace_all(dset_name, "\\.", "_"), ".rds"))

