
# Subplot -----------------------------------------------------------------

library(tidyverse)
library(stringr)
library(broom)
library(survival)

source("utilities-downstream.R")

if (str_sub(dset_name, -2, -1) %in% c("sl", "sn")) {
    str_censoredInt <- "0"
    cen_val <- 0
} else {
    str_censoredInt <- "NA"
    cen_val <- NA
}

load(paste0("fltr/", dset_name, ".rda"))
eval(parse(text = paste0("df_mss <- ", dset_name)))
df_orig <- modify_badftrpk(df_mss)

# Filtering
res_orig <- gen_subplot(df_orig, censoredInt = str_censoredInt)

# Subplot summarization
saveRDS(res_orig, paste0("subplot/", str_replace_all(dset_name, "\\.", "_"), ".rds"))

vers <- c("v2", "v3", "v4", "v5", "v6", "v7")

for (i in seq_along(vers)) {
    meth_ver <- vers[i]
    load(paste0("fltr/", meth_ver, "/", dset_name, ".rmftr.rda"))
    load(paste0("fltr/", meth_ver, "/", dset_name, ".rmpk.rda"))
    eval(parse(text = paste0("df_rmftr <- ", dset_name, ".rmftr")))
    eval(parse(text = paste0("df_rmpk <- ", dset_name, ".rmpk")))
    df_fltr <- modify_badftrpk(df_mss, df_rmftr, df_rmpk, cen_val = cen_val)
    
    res_fltr <- gen_subplot(df_fltr, censoredInt = str_censoredInt)
    saveRDS(res_fltr, paste0("subplot/", str_replace_all(dset_name, "\\.", "_"), "_m5", meth_ver, ".rds"))
}

rm(list = ls())
