
# Load --------------------------------------------------------------------

library(tidyverse)
library(stringr)

# dset_name <- "dia.bruderer.sn"
dset_name <- "dia.bruderer.sl"

meth <- c("All", "Top10", "Top5", "Top3", "Method 1", "Method 4")
meth_code <- c("all", "top10", "top5", "top3", "m1.v1", "m4.v2")

test_res <- vector("list", length(meth_code))
for (i in seq_along(meth_code)) {
    test_name <- paste0("test.", dset_name, ".", meth_code[i])
    load(paste0("test/", test_name, ".rda"))
    eval(parse(text = paste0("test_res[[i]] <- ", test_name)))
    test_res[[i]] <- test_res[[i]] %>% as_tibble() %>% mutate(method = meth[i])
}
test_res <- bind_rows(test_res)


# Ground truth ------------------------------------------------------------

mix <- list(
    c("P02754", "P00921", "P80025", "P02662", "P00366"), 
    c("P12799", "P02672", "P02789", "P02676", "P61823"), 
    c("P68082", "P02666")
)

conc <- list(
    c(1.5, 1.65, 1.815, 1.995, 15, 16.515, 18.165, 19.995), 
    c(100, 62.995, 39.685, 25, 2, 1.26, 0.795, 0.5),
    c(0.05, 0.2, 0.8, 3.2, 12.8, 51.2, 204.8, 819.2)
)

cond_code <- "S"
nb_cond <- length(conc[[1]])
lab <- rep(NA, (nb_cond ^ 2 - nb_cond) / 2)
fc <- vector("list", length(conc))
for (m in 1:length(fc)) {
    fc[[m]] <- rep(NA, (nb_cond ^ 2 - nb_cond) / 2)
}

idx <- 1
for (i in 1:(nb_cond - 1)) {
    for (j in (i + 1):nb_cond) {
        lab[idx] <- paste0(cond_code, j, "-", cond_code, i)
        for (m in 1:length(fc)) {
            fc[[m]][idx] <- conc[[m]][j] / conc[[m]][i]
        }
        idx <- idx + 1
    }
}

gold <- tibble(
    Protein = rep(unlist(mix), each = length(lab)), 
    Label = rep(lab, sum(map_int(mix, length))), 
    trueFC = unlist(map2(fc, map(mix, length), rep))
)

# mx1 <- c("P02754", "P00921", "P80025", "P02662", "P00366")
# mx2 <- c("P12799", "P02672", "P02789", "P02676", "P61823")
# mx3 <- c("P68082", "P02666")
# 
# conc1 <- c(1.5, 1.65, 1.815, 1.995, 15, 16.515, 18.165, 19.995)
# conc2 <- c(100, 62.995, 39.685, 25, 2, 1.26, 0.795, 0.5)
# conc3 <- c(0.05, 0.2, 0.8, 3.2, 12.8, 51.2, 204.8, 819.2)
# 
# nb_cond <- 8
# lab <- fc1 <- fc2 <- fc3 <- rep(NA, (nb_cond ^ 2 - nb_cond) / 2)
# idx <- 1
# for (i in 1:(nb_cond - 1)) {
#     for (j in (i + 1):nb_cond) {
#         lab[idx] <- paste0("S", j, "-S", i)
#         fc1[idx] <- conc1[j] / conc1[i]
#         fc2[idx] <- conc2[j] / conc2[i]
#         fc3[idx] <- conc3[j] / conc3[i]
#         idx <- idx + 1
#     }
# }
# 
# gold <- bind_rows(
#     tibble(Protein = rep(mx1, each = length(lab)), Label = rep(lab, length(mx1)), trueFC = rep(fc1, length(mx1))), 
#     tibble(Protein = rep(mx2, each = length(lab)), Label = rep(lab, length(mx2)), trueFC = rep(fc2, length(mx2))), 
#     tibble(Protein = rep(mx3, each = length(lab)), Label = rep(lab, length(mx3)), trueFC = rep(fc3, length(mx3)))
# )

test_res <- bind_rows(
    test_res %>% filter(Protein %in% unlist(mix)) %>% left_join(gold), 
    test_res %>% filter(!(Protein %in% unlist(mix))) %>% mutate(trueFC = 1)
)

test_res <- test_res %>% 
    mutate(de = ifelse(trueFC == 1, "No change", "With change")) %>% 
    mutate(trueLog2FC = log2(trueFC), squrError = (log2FC - trueLog2FC) ^ 2)

only_fc2 <- FALSE
if (only_fc2) {
    test_res <- test_res %>% 
        filter(abs(trueLog2FC) <= 1)
}


# Summarize estimation & testing ------------------------------------------

sum_estimation <- test_res %>% 
    filter(!is.na(pvalue)) %>%
    group_by(method, de) %>% 
    summarise(ave_se = mean(squrError), med_se = median(squrError)) %>% 
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


sensPrecCurves <- function(test_results) {
    meths <- unique(test_results$method)
    test_meth = vector("list", length = length(meths))
    for (i in seq_along(meths)) {
        test_sub <- test_results %>% 
            filter(method == meths[i], !is.na(pvalue)) %>% 
            arrange(adj.pvalue)
        ordered_de <- test_sub$de == "With change"
        test_meth[[i]] <- tibble(
            sens = cumsum(ordered_de) / sum(ordered_de), 
            fdr = cumsum(1 - ordered_de) / seq_along(ordered_de),
            fpr = cumsum(!ordered_de) / sum(!ordered_de)
        ) %>% mutate(method = meths[i])
    }
    
    return(bind_rows(test_meth))
}

all_sp_curves <- sensPrecCurves(test_res)


# Examine results ---------------------------------------------------------

# Precision-recall
all_sp_curves %>% 
    mutate(method = factor(method, levels = meth)) %>% 
    ggplot(aes(x = fdr, y = sens, group = method, color = method)) + 
    geom_line() + 
    coord_cartesian(xlim = c(0, 0.15), ylim = c(0.25, 1)) + 
    ggtitle(paste0("Sensitivity vs. FDR (", dset_name, ")"))
ggsave(paste0("prc-", str_replace_all(dset_name, "\\.", "_"), ".png"), width = 6, height = 4)

# ROC curve
all_sp_curves %>% 
    mutate(method = factor(method, levels = meth)) %>% 
    ggplot(aes(x = fpr, y = sens, group = method, color = method)) + 
    geom_line() + 
    coord_cartesian(xlim = c(0, 0.15), ylim = c(0.4, 1)) + 
    ggtitle(paste0("ROC curve (", dset_name, ")"))
ggsave(paste0("roc-", str_replace_all(dset_name, "\\.", "_"), ".png"), width = 6, height = 4)

sum_testing

sum_estimation %>% select(-ave_se) %>% spread(de, med_se)
sum_estimation %>% select(-med_se) %>% spread(de, ave_se)

test_res %>%
    filter(abs(log2FC) != Inf) %>%
    ggplot(aes(x = method, y = sqrt(squrError))) +
    geom_boxplot() +
    facet_wrap(~ de) +
    coord_cartesian(ylim = c(0, 1.5))

test_res %>%
    filter(abs(log2FC) != Inf) %>% filter(Protein %in% mx1) %>% filter(trueLog2FC < 2) %>% 
    ggplot(aes(x = trueLog2FC, y = sqrt(squrError), color = method, group = method)) +
    geom_point(position = position_dodge(width = 0.075), alpha = 0.5, size = 2)

test_res %>%
    filter(abs(log2FC) != Inf) %>% filter(Protein %in% mx1) %>% filter(trueLog2FC > 2) %>% 
    ggplot(aes(x = trueLog2FC, y = sqrt(squrError), color = method, group = method)) +
    geom_point(position = position_dodge(width = 0.075), alpha = 0.5, size = 2)

test_res %>% 
    filter(abs(log2FC) != Inf) %>% 
    ggplot(aes(method, SE, color = method)) + 
    geom_boxplot() + 
    facet_wrap(~ de)

