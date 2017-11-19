
# Cross-tool comparison ---------------------------------------------------

library(tidyverse)
library(stringr)
library(broom)

source("utilities-downstream.R")

# expt_name <- "dia.bruderer"
# expt_name <- "dia.navarro"
# expt_name <- "dia.selevsek"
# 
# expt_name <- "dda.choi"
# expt_name <- "dda.iprg"

alldset_name <- str_replace(dir("fltr", pattern = expt_name), ".rda", "")

meth_ver <- "m5v2"

list_orig <- vector("list", length(alldset_name))
list_fltr <- vector("list", length(alldset_name))
list_orig_sub <- vector("list", length(alldset_name))
list_fltr_sub <- vector("list", length(alldset_name))

regex_uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")

for (i in seq_along(alldset_name)) {
    load(paste0("fltr/v2/", alldset_name[i], ".rmftr.rda"))
    load(paste0("fltr/v2/", alldset_name[i], ".rmpk.rda"))
    eval(parse(text = paste0("df_rmftr <- ", alldset_name[i], ".rmftr")))
    eval(parse(text = paste0("df_rmpk <- ", alldset_name[i], ".rmpk")))
    prot_mod <- union(unique(df_rmftr$protein), unique(df_rmpk$protein))
    # Subplot summarization
    dset_name_us <- str_replace_all(alldset_name[i], "\\.", "_")
    list_orig[[i]] <- readRDS(paste0("subplot/", dset_name_us, ".rds"))$res_subplot %>% 
        mutate(sw = str_to_upper(str_sub(dset_name_us, -2, -1)))
    list_fltr[[i]] <- readRDS(paste0("subplot/", dset_name_us, "_", meth_ver, ".rds"))$res_subplot %>% 
        mutate(sw = str_to_upper(str_sub(dset_name_us, -2, -1)))
    # Keep only those altered proteins 
    list_orig_sub[[i]] <- list_orig[[i]] %>% filter(protein %in% prot_mod)
    list_fltr_sub[[i]] <- list_fltr[[i]] %>% filter(protein %in% prot_mod)
}

if (expt_name == "dia.selevsek") {
    df_orig <- bind_rows(list_orig)
    df_fltr <- bind_rows(list_fltr)
    df_orig_sub <- bind_rows(list_orig_sub)
    df_fltr_sub <- bind_rows(list_fltr_sub)
} else {
    df_orig <- bind_rows(list_orig) %>% 
        mutate(protein = str_extract(protein, pattern = regex_uniprot_iso))
    df_fltr <- bind_rows(list_fltr) %>% 
        mutate(protein = str_extract(protein, pattern = regex_uniprot_iso))
    df_orig_sub <- bind_rows(list_orig_sub) %>% 
        mutate(protein = str_extract(protein, pattern = regex_uniprot_iso))
    df_fltr_sub <- bind_rows(list_fltr_sub) %>% 
        mutate(protein = str_extract(protein, pattern = regex_uniprot_iso))
}


df_orig %>% spread(sw, log2inty) %>% select(-protein, -run) %>% 
    cor(., use = "pairwise.complete.obs")

df_fltr %>% spread(sw, log2inty) %>% select(-protein, -run) %>% 
    cor(., use = "pairwise.complete.obs")

df_orig_sub %>% spread(sw, log2inty) %>% select(-protein, -run) %>% 
    cor(., use = "pairwise.complete.obs")

df_fltr_sub %>% spread(sw, log2inty) %>% select(-protein, -run) %>% 
    cor(., use = "pairwise.complete.obs")

