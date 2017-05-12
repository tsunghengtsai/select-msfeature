# Load required packages --------------------------------------------------
library(MASS)
library(tidyverse)
library(stringr)
library(broom)


# Simulation data ---------------------------------------------------------
nb_sim <- 500
sizes <- 2:100

df_sim <- tibble(
    size = rep(sizes, each = nb_sim), 
    sim = rep(1:nb_sim, length(sizes))
)

df_sim <- df_sim %>% 
    mutate(data = map(size, rnorm))

df_sim <- df_sim %>% 
    mutate(
        samp_sd = map_dbl(data, sd),
        samp_mad = map_dbl(data, mad)
    )

# Summary of the simulation
df_sum <- df_sim %>% 
    group_by(size) %>% 
    summarise(
        ave_sd = mean(samp_sd), 
        std_sd = sd(samp_sd), 
        ave_mad = mean(samp_mad), 
        std_mad = sd(samp_mad)
    )

ave_sum <- df_sum %>% 
    select(-std_sd, -std_mad) %>% 
    gather(estimator, ave_sim, ave_sd, ave_mad) %>% 
    mutate(estimator = ifelse(str_detect(estimator, "sd"), "SD", "MAD"))

sd_sum <- df_sum %>% 
    select(-ave_sd, -ave_mad) %>% 
    gather(estimator, sd_sim, std_sd, std_mad) %>% 
    mutate(estimator = ifelse(str_detect(estimator, "sd"), "SD", "MAD"))

df_sum <- left_join(ave_sum, sd_sum)

# Mean +/- SD of the estimates
df_sum %>% 
    ggplot(aes(size, ave_sim, ymin = ave_sim - sd_sim, ymax = ave_sim + sd_sim, 
               colour = estimator, fill = estimator)) + 
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    geom_hline(yintercept = 1) + 
    labs(y = "estimate", title = str_c(nb_sim, " simulations"))

# SD of sample estimates
df_sum %>% 
    ggplot(aes(size, sd_sim, colour = estimator, fill = estimator)) + 
    geom_point(position = position_dodge(width = 0.7)) + 
    geom_hline(yintercept = filter(df_sum, estimator == "MAD", size %in% c(10, 50, 100)) %>% .$sd_sim, lty = 2) + 
    labs(y = "estimate", title = str_c(nb_sim, " simulations"))


# Functions for two-way analysis ------------------------------------------

gen2way <- function(n, p) {
    df_2way <- tibble(
        run = rep(str_pad(1:n, str_length(n), pad = "0"), each = p), 
        ftr = rep(str_pad(1:p, str_length(p), pad = "0"), n), 
        obs = rnorm(n * p)
    )
    
    return(df_2way)
}

gen2way_olr <- function(n, p) {
    df_2way <- tibble(
        run = rep(str_pad(1:n, str_length(n), pad = "0"), each = p), 
        ftr = rep(str_pad(1:p, str_length(p), pad = "0"), n), 
        obs = c(20, rnorm(n * p - 1))
    )
    
    return(df_2way)
}

mp_df <- function(df) {
    obs_wide <- df %>% spread(ftr, obs)
    obs_mat <- data.matrix(obs_wide[, -1])
    rownames(obs_mat) <- obs_wide$run
    mp_fit <- medpolish(obs_mat, na.rm = TRUE, trace.iter = FALSE)
    
    return(mp_fit)
}

mad_mpres <- function(mp_fit) {
    res <- mp_fit$residuals %>% as.vector()

    return(mad(res))
}

adjmad_mpres <- function(mp_fit, dfree) {
    res <- mp_fit$residuals %>% as.vector()
    adjmad <- mad(res) * sqrt((length(res) - 1) / dfree)
    
    return(adjmad)
}


# Two-way analysis --------------------------------------------------------

# Simulation setting
nb_sim <- 50
nb_run <- 2 ^ (1:5)
nb_ftr <- c(3, 10, 20)


# Fitting with lm and rlm -------------------------------------------------
library(robustbase)
library(L1pack)

df_sim2 <- tibble(
    n_run = rep(nb_run, each = length(nb_ftr) * nb_sim), 
    n_ftr = rep(nb_ftr, length(nb_run) * nb_sim), 
    sim = rep(rep(1:nb_sim, each = length(nb_ftr)), length(nb_run))
)

df_sim2 <- df_sim2 %>% 
    mutate(data = map2(n_run, n_ftr, gen2way))

df_sim2 <- df_sim2 %>% 
    mutate(data = map2(n_run, n_ftr, gen2way_olr))


df_sim2 <- df_sim2 %>%
    mutate(
        lm_fit = map(data, ~ lm(obs ~ run + ftr, data = .)),
        rlm_fit = map(data, ~ rlm(obs ~ run + ftr, data = .)), 
        mp_fit = map(data, mp_df), 
        lad_fit = map(data, ~ lad(obs ~ run + ftr, data = ., method = "EM")),
        rob_fit = map(data, ~ lmrob(obs ~ run + ftr, data = .))
    ) %>% 
    mutate(
        samp_lms = map_dbl(lm_fit, sigma), 
        samp_rlms = map_dbl(rlm_fit, ~ summary(.)$sigma), 
        dfree = map_int(lm_fit, df.residual), 
        samp_mad = map_dbl(mp_fit, mad_mpres),
        samp_mps = map2_dbl(mp_fit, dfree, adjmad_mpres), 
        samp_lads = map_dbl(lad_fit, ~ .$scale), 
        samp_robs = map_dbl(rob_fit, ~ summary(.)$sigma)
    )

# Summary of the simulation
df_sum2 <- df_sim2 %>% 
    group_by(n_ftr, n_run) %>% 
    summarise(
        ave_lms = mean(samp_lms), 
        std_lms = sd(samp_lms), 
        ave_rlms = mean(samp_rlms), 
        std_rlms = sd(samp_rlms), 
        ave_mps = mean(samp_mps), 
        std_mps = sd(samp_mps), 
        ave_lads = mean(samp_lads), 
        std_lads = sd(samp_lads), 
        ave_robs = mean(samp_robs),
        std_robs = sd(samp_robs)
    ) %>% 
    ungroup()

ave_sum2 <- df_sum2 %>% 
    dplyr::select(-std_lms, -std_rlms, -std_mps, -std_lads, -std_robs) %>% 
    gather(estimator, ave_sim, ave_lms, ave_rlms, ave_mps, ave_lads, ave_robs) %>% 
    mutate(estimator = str_replace(estimator, "ave_", ""))

sd_sum2 <- df_sum2 %>% 
    dplyr::select(-ave_lms, -ave_rlms, -ave_mps, -ave_lads, -ave_robs) %>% 
    gather(estimator, sd_sim, std_lms, std_rlms, std_mps, std_lads, std_robs) %>% 
    mutate(estimator = str_replace(estimator, "std_", ""))

df_sum2_td <- left_join(ave_sum2, sd_sum2)

# Mean +/- SD of the estimates
df_sum2_td %>% 
    ggplot(aes(log2(n_run), ave_sim, ymin = ave_sim - sd_sim, ymax = ave_sim + sd_sim, 
               colour = estimator, fill = estimator)) + 
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    geom_hline(yintercept = 1) + 
    coord_cartesian(ylim = c(0, 2)) + 
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) + 
    facet_wrap(~ n_ftr)

# SD of sample estimates
df_sum2_td %>% 
    ggplot(aes(log2(n_run), sd_sim, colour = estimator, fill = estimator)) + 
    geom_point(position = position_dodge(width = 0.7)) + 
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) + 
    facet_wrap(~ n_ftr)



# Compare scale estimate in rlm -------------------------------------------

nb_sim <- 100
nb_run <- 2 ^ (1:6)
nb_ftr <- c(3, 10, 20, 30)

df_sim3 <- tibble(
    n_run = rep(nb_run, each = length(nb_ftr) * nb_sim), 
    n_ftr = rep(nb_ftr, length(nb_run) * nb_sim), 
    sim = rep(rep(1:nb_sim, each = length(nb_ftr)), length(nb_run))
)

# df_sim3 <- df_sim3 %>% 
#     mutate(data = map2(n_run, n_ftr, gen2way))

df_sim3 <- df_sim3 %>%
    mutate(data = map2(n_run, n_ftr, gen2way_olr))


df_sim3 <- df_sim3 %>%
    mutate(
        lm_fit = map(data, ~ lm(obs ~ run + ftr, data = .)),
        rlm_fit0 = map(data, ~ rlm(obs ~ run + ftr, data = .)), 
        rlm_fit1 = map(data, ~ rlm(obs ~ run + ftr, data = ., scale.est = "Huber")),
        rlm_fit2 = map(data, ~ rlm(obs ~ run + ftr, data = ., scale.est = "proposal 2"))
    ) %>% 
    mutate(
        samp_lms = map_dbl(lm_fit, sigma), 
        samp_rlms0 = map_dbl(rlm_fit0, ~ summary(.)$sigma), 
        samp_rlms1 = map_dbl(rlm_fit1, ~ summary(.)$sigma), 
        samp_rlms2 = map_dbl(rlm_fit2, ~ summary(.)$sigma)
    )


# Summary of the simulation
df_sum3 <- df_sim3 %>% 
    group_by(n_ftr, n_run) %>% 
    summarise(
        ave_lms = mean(samp_lms), 
        std_lms = sd(samp_lms), 
        ave_rlms0 = mean(samp_rlms0), 
        std_rlms0 = sd(samp_rlms0), 
        ave_rlms1 = mean(samp_rlms1), 
        std_rlms1 = sd(samp_rlms1), 
        ave_rlms2 = mean(samp_rlms2), 
        std_rlms2 = sd(samp_rlms2)
    ) %>% 
    ungroup()

ave_sum3 <- df_sum3 %>% 
    dplyr::select(-std_lms, -std_rlms0, -std_rlms1, -std_rlms2) %>% 
    gather(estimator, ave_sim, ave_lms, ave_rlms0, ave_rlms1, ave_rlms2) %>% 
    mutate(estimator = str_replace(estimator, "ave_", ""))

sd_sum3 <- df_sum3 %>% 
    dplyr::select(-ave_lms, -ave_rlms0, -ave_rlms1, -ave_rlms2) %>% 
    gather(estimator, sd_sim, std_lms, std_rlms0, std_rlms1, std_rlms2) %>% 
    mutate(estimator = str_replace(estimator, "std_", ""))

df_sum3_td <- left_join(ave_sum3, sd_sum3)

# Mean +/- SD of the estimates
df_sum3_td %>% 
    ggplot(aes(log2(n_run), ave_sim, ymin = ave_sim - sd_sim, ymax = ave_sim + sd_sim, 
               colour = estimator, fill = estimator)) + 
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    geom_hline(yintercept = 1) + 
    coord_cartesian(ylim = c(0, 2)) + 
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) + 
    facet_wrap(~ n_ftr)

# SD of sample estimates
df_sum3_td %>%
    ggplot(aes(log2(n_run), sd_sim, colour = estimator, fill = estimator)) +
    geom_point(position = position_dodge(width = 0.7)) +
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) +
    facet_wrap(~ n_ftr)
