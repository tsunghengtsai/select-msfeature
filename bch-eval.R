
library(tidyverse)
library(stringr)
library(broom)

source("utilities-downstream.R")

test_res <- readRDS(paste0("esttest/", str_replace_all(dset_name, "\\.", "_"), ".rds"))

# Gold standard
if (str_detect(dset_name, "choi")) {
    source("gold_dda_choi.R")
    
    test_res <- bind_rows(
        test_res %>% 
            filter(protein %in% unlist(mix)) %>% 
            left_join(gold %>% rename(protein = Protein, ctrx = Label)), 
        test_res %>% filter(!(protein %in% unlist(mix))) %>% mutate(trueFC = 1)
    ) %>% 
        mutate(de = ifelse(trueFC == 1, "No change", "With change")) %>% 
        mutate(trueLog2FC = log2(trueFC), est_err2 = (logFC - trueLog2FC) ^ 2)
    
} else if (str_detect(dset_name, "iprg")) {
    source("gold_dda_iprg.R")
    
    test_res <- bind_rows(
        test_res %>% 
            filter(protein %in% unlist(mix)) %>% 
            left_join(gold %>% rename(protein = Protein, ctrx = Label)), 
        test_res %>% filter(!(protein %in% unlist(mix))) %>% mutate(trueFC = 1)
    ) %>% 
        mutate(de = ifelse(trueFC == 1, "No change", "With change")) %>% 
        mutate(trueLog2FC = log2(trueFC), est_err2 = (logFC - trueLog2FC) ^ 2)
    
} else if (str_detect(dset_name, "navarro")) {
    
    gold <- test_res %>% 
        distinct(protein) %>% 
        mutate(trueFC = ifelse(str_detect(protein, "HUMAN"), 1, ifelse(str_detect(protein, "ECOLI"), 0.25, 2)))
    
    test_res <- test_res %>% 
        left_join(gold) %>% 
        mutate(de = ifelse(trueFC == 1, "No change", "With change")) %>% 
        mutate(trueLog2FC = log2(trueFC), est_err2 = (logFC - trueLog2FC) ^ 2)

}


all_sp_curves <- sensPrecCurves(test_res)

sum_estimation <- test_res %>% 
    filter(!is.na(pvalue)) %>%
    group_by(method, de) %>% 
    summarise(ave_se = mean(est_err2), med_se = median(est_err2)) %>% 
    ungroup()

sum_testing <- test_res %>% 
    filter(!is.na(pvalue)) %>%
    group_by(method) %>% 
    summarise(
        n_test = n(), n_pos = sum(de == "With change"), n_neg = sum(de == "No change"), 
        n_ppos = sum(adj.pvalue <= 0.05), n_pneg = sum(adj.pvalue > 0.05),
        n_tp = sum(adj.pvalue <= 0.05 & de == "With change"), 
        n_tn = sum(adj.pvalue > 0.05 & de == "No change"), 
        n_fp = sum(adj.pvalue <= 0.05 & de == "No change")
    ) %>% 
    mutate(
        sens = n_tp / n_pos, 
        spec = n_tn / n_neg, 
        ppv = n_tp / (n_tp + n_fp), 
        fdr = n_fp / (n_tp + n_fp)
    )

all_sp_curves %>% 
    # mutate(method = factor(method, levels = meth)) %>% 
    ggplot(aes(x = fdr, y = sens, group = method, color = method)) + 
    geom_line() + 
    coord_cartesian(xlim = c(0, 0.15), ylim = c(0.25, 1)) + 
    ggtitle(paste0("Sensitivity vs. FDR (", dset_name, ")"))

all_sp_curves %>% 
    # mutate(method = factor(method, levels = meth)) %>% 
    ggplot(aes(x = fpr, y = sens, group = method, color = method)) + 
    geom_line() + 
    coord_cartesian(xlim = c(0, 0.15), ylim = c(0.7, 1)) + 
    ggtitle(paste0("ROC curve (", dset_name, ")"))

sum_testing
sum_estimation %>% select(method, de, ave_se) %>% spread(de, ave_se)
sum_estimation %>% select(method, de, med_se) %>% spread(de, med_se)
