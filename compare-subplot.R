
# Load data and carry out subplot summarization ---------------------------

library(tidyverse)
library(stringr)
library(survival)

source("utilities.R")

dset_name <- "dda.choi.mq"
# dset_name <- "dda.choi.pd"
# dset_name <- "dda.choi.pg"
# dset_name <- "dda.choi.sl"
# dset_name <- "dia.navarro.os"
# dset_name <- "dia.navarro.sl"
# dset_name <- "dia.navarro.sn"

# meth_code <- "m5_v1"
meth_code <- "m5_v2"
# meth_code <- "m5_v3"

# input
load(paste0("case/", dset_name, ".rda"))
load(paste0("case/", meth_code, "/", dset_name, ".rmftr.rda"))
load(paste0("case/", meth_code, "/", dset_name, ".rmpk.rda"))

eval(parse(text = paste0("df_all <- ", dset_name)))
eval(parse(text = paste0("df_rmftr <- ", dset_name, ".rmftr")))
eval(parse(text = paste0("df_rmpk <- ", dset_name, ".rmpk")))

rm(list = paste0(dset_name, c("", ".rmftr", ".rmpk")))

if (str_sub(dset_name, -2, -1) %in% c("sl", "sn")) {
    str_censoredInt <- "0"
    cen_val <- 0
} else {
    str_censoredInt <- "NA"
    cen_val <- NA
}

case_name <- paste0("checklist_se_top10_", str_replace_all(dset_name, "\\.", "_"))
load(paste0("case/", case_name, ".rda"))
case_prots <- unique(top.issue.protein$Protein)

df_orig <- modify_badftrpk(df_all, cen_val = cen_val)
df_fltr <- modify_badftrpk(df_all, df_rmftr, df_rmpk, cen_val = cen_val)
df_rmolr <- modify_badftrpk(df_all, df_rmftr, df_rmpk, cen_val = cen_val, rm_olr = TRUE)

res_orig <- gen_subplot(df_orig %>% filter(protein %in% case_prots), censoredInt = str_censoredInt)
res_fltr <- gen_subplot(df_fltr %>% filter(protein %in% case_prots), censoredInt = str_censoredInt)
res_rmolr <- gen_subplot(df_rmolr %>% filter(protein %in% case_prots), censoredInt = str_censoredInt)

subplot <- bind_rows(
    res_orig$res_subplot %>% mutate(meth = "All"), 
    res_fltr$res_subplot %>% mutate(meth = "Filtered"), 
    res_rmolr$res_subplot %>% mutate(meth = "Removed")
) %>% 
    left_join(df_orig %>% distinct(run, subject, group)) %>% 
    mutate(run = str_pad(run, 2, side = "left", pad = "0")) %>% 
    mutate(grp_run = str_c(group, run, sep = "_"))

aft <- bind_rows(
    res_orig$res_aft %>% mutate(meth = "All"), 
    res_fltr$res_aft %>% mutate(meth = "Filtered"), 
    res_rmolr$res_aft %>% mutate(meth = "Removed")
) %>% 
    select(protein, feature:subject, ind_obs, meth) %>% 
    mutate(run = str_pad(run, 2, side = "left", pad = "0")) %>% 
    mutate(grp_run = str_c(group, run, sep = "_"))



# Plots for comparison ----------------------------------------------------

plot_aftprofile <- function(aftdata, protein_name) {
    oneaft <- aft %>% filter(protein == protein_name)
    oneaft %>% 
        ggplot(aes(run, log2inty, group = feature, color = feature)) + 
        geom_point(aes(alpha = factor(!is_censored)), size = 2) + 
        scale_alpha_discrete(range = c(0.4, 1)) + 
        geom_line() + 
        facet_wrap(~ meth) + 
        labs(title = "Profile plot after AFT imputation", subtitle = protein_name) + 
        guides(alpha = FALSE, color = FALSE)
}
plot_subaft <- function(aftdata, subdata, protein_name) {
    oneaft <- aft %>% filter(protein == protein_name)
    onesub <- subplot %>% filter(protein == protein_name)
    oneaft %>% 
        ggplot(aes(run, log2inty, group = feature)) + 
        geom_point(aes(alpha = factor(!is_censored)), size = 2, color = "gray") + 
        scale_alpha_discrete(range = c(0.4, 1)) + 
        geom_line(color = "gray") + 
        geom_point(data = onesub, aes(run, log2inty, group = meth, color = meth)) + 
        geom_line(data = onesub, aes(run, log2inty, group = meth, color = meth)) + 
        facet_wrap(~ meth) + 
        labs(title = "Subplot summarization", subtitle = protein_name) + 
        guides(alpha = FALSE, color = FALSE)
}

pdf(paste0("case_", str_replace_all(dset_name, "\\.", "_"), ".pdf"), width = 12, height = 6)
for (i in seq_along(case_prots)) {
    print(plot_aftprofile(aft, case_prots[i]))
    print(plot_subaft(aft, subplot, case_prots[i]))
}
dev.off()


