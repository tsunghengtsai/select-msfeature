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


# Two-way analysis with different estimators ------------------------------
library(L1pack)

nb_sim <- 100
nb_run <- 2 ^ (1:6)
nb_ftr <- c(3, 5, 10, 20)

df_sim2 <- tibble(
    n_run = rep(nb_run, each = length(nb_ftr) * nb_sim), 
    n_ftr = rep(nb_ftr, length(nb_run) * nb_sim), 
    sim = rep(rep(1:nb_sim, each = length(nb_ftr)), length(nb_run))
)

df_sim2 <- df_sim2 %>%
    mutate(data = map2(n_run, n_ftr, gen2way))

# df_sim2 <- df_sim2 %>%
#     mutate(data = map2(n_run, n_ftr, gen2way_olr))


df_sim2 <- df_sim2 %>%
    mutate(
        lm_fit = map(data, ~ lm(obs ~ run + ftr, data = .)),
        rlm_fit = map(data, ~ rlm(obs ~ run + ftr, data = .)), 
        rlm_fit_huber = map(data, ~ rlm(obs ~ run + ftr, data = ., scale.est = "Huber")),
        mp_fit = map(data, mp_df), 
        lad_fit = map(data, ~ lad(obs ~ run + ftr, data = ., method = "EM"))
    ) %>% 
    mutate(
        s_lm = map_dbl(lm_fit, sigma), 
        s_rlm_mad = map_dbl(rlm_fit, ~ summary(.)$sigma), 
        s_rlm_huber = map_dbl(rlm_fit_huber, ~ summary(.)$sigma), 
        s_mp = map_dbl(mp_fit, mad_mpres),
        s_lad = map_dbl(lad_fit, ~ .$scale), 
        s_lad_mad = map_dbl(lad_fit, ~ mad(.$residuals))
    )


# Summary of the simulation
df_sum2 <- df_sim2 %>% 
    group_by(n_ftr, n_run) %>% 
    summarise_each(funs(mean, sd), s_lm:s_lad_mad) %>% 
    ungroup() %>% 
    gather(est, value, s_lm_mean:s_lad_mad_sd) %>% 
    mutate(summary = ifelse(str_detect(est, "_mean"), "sim_mean", "sim_sd")) %>% 
    mutate(est = str_replace_all(est, "s_|_mean|_sd", "")) %>% 
    spread(summary, value)
    

# Mean +/- SD of the estimates
df_sum2 %>% 
    ggplot(aes(as.factor(n_run), sim_mean, ymin = sim_mean - sim_sd, ymax = sim_mean + sim_sd, colour = est, fill = est)) + 
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    geom_hline(yintercept = 1) + 
    coord_cartesian(ylim = c(0, 2)) +
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) + 
    facet_wrap(~ n_ftr)


# Compare scale estimate in rlm -------------------------------------------

nb_sim <- 100
nb_run <- 2 ^ (1:6)
nb_ftr <- c(3, 5, 10, 20)

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
        s_lm = map_dbl(lm_fit, sigma), 
        s_rlm0 = map_dbl(rlm_fit0, ~ summary(.)$sigma), 
        s_rlm1 = map_dbl(rlm_fit1, ~ summary(.)$sigma), 
        s_rlm2 = map_dbl(rlm_fit2, ~ summary(.)$sigma)
    )


# Summary of the simulation
df_sum3 <- df_sim3 %>% 
    group_by(n_ftr, n_run) %>% 
    summarise_each(funs(mean, sd), s_lm:s_rlm2) %>% 
    ungroup() %>% 
    gather(est, value, s_lm_mean:s_lad_mad_sd) %>% 
    mutate(summary = ifelse(str_detect(est, "_mean"), "sim_mean", "sim_sd")) %>% 
    mutate(est = str_replace_all(est, "s_|_mean|_sd", "")) %>% 
    spread(summary, value)


# Mean +/- SD of the estimates
df_sum3 %>% 
    ggplot(aes(as.factor(n_run), sim_mean, ymin = sim_mean - sim_sd, ymax = sim_mean + sim_sd, colour = est, fill = est)) + 
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    geom_hline(yintercept = 1) + 
    coord_cartesian(ylim = c(0, 2)) +
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) + 
    facet_wrap(~ n_ftr)

# SD of sample estimates
df_sum3 %>%
    ggplot(aes(as.factor(n_run), sim_sd, colour = est, fill = est)) +
    geom_point(position = position_dodge(width = 0.7)) +
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) +
    facet_wrap(~ n_ftr)

