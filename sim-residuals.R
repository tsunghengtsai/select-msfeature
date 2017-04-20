# Load required packages --------------------------------------------------
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


# Data in two-way table ---------------------------------------------------
gen2way <- function(n, p) {
    df_2way <- tibble(
        run = rep(1:n, each = p), 
        ftr = rep(1:p, n), 
        obs = rnorm(n * p)
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


# Simulation setting
nb_sim <- 100
nb_run <- 2:100
nb_ftr <- c(3, 10, 20)

df_sim2 <- tibble(
    n_run = rep(nb_run, each = length(nb_ftr) * nb_sim), 
    n_ftr = rep(nb_ftr, length(nb_run) * nb_sim), 
    sim = rep(rep(1:nb_sim, each = length(nb_ftr)), length(nb_run))
)

# Generate simulation data
df_sim2 <- df_sim2 %>% 
    mutate(data = map2(n_run, n_ftr, gen2way))

# Fit lm and medpolish models
df_sim2 <- df_sim2 %>% 
    mutate(
        lm_fit = map(data, ~ lm(obs ~ run + ftr, data = .)), 
        mp_fit = map(data, mp_df)
    )

# Estimates with lm and medpolish
df_sim2 <- df_sim2 %>% 
    mutate(
        samp_lms = map_dbl(lm_fit, sigma), 
        dfree = map_int(lm_fit, df.residual), 
        samp_mps = map_dbl(mp_fit, mad_mpres)
    )

# Summary of the simulation
df_sum2 <- df_sim2 %>% 
    group_by(n_ftr, n_run) %>% 
    summarise(
        ave_lms = mean(samp_lms), 
        std_lms = sd(samp_lms), 
        ave_mps = mean(samp_mps), 
        std_mps = sd(samp_mps)
    )

ave_sum2 <- df_sum2 %>% 
    select(-std_lms, -std_mps) %>% 
    gather(estimator, ave_sim, ave_lms, ave_mps) %>% 
    mutate(estimator = ifelse(str_detect(estimator, "lms"), "LM", "MP-MAD"))

sd_sum2 <- df_sum2 %>% 
    select(-ave_lms, -ave_mps) %>% 
    gather(estimator, sd_sim, std_lms, std_mps) %>% 
    mutate(estimator = ifelse(str_detect(estimator, "lms"), "LM", "MP-MAD"))

df_sum2 <- left_join(ave_sum2, sd_sum2)

# Mean +/- SD of the estimates
df_sum2 %>% 
    ggplot(aes(n_run, ave_sim, ymin = ave_sim - sd_sim, ymax = ave_sim + sd_sim, 
               colour = estimator, fill = estimator)) + 
    geom_pointrange(position = position_dodge(width = 0.7)) + 
    geom_hline(yintercept = 1) + 
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) + 
    facet_wrap(~ n_ftr)

# SD of sample estimates
df_sum2 %>% 
    ggplot(aes(n_run, sd_sim, colour = estimator, fill = estimator)) + 
    geom_point(position = position_dodge(width = 0.7)) + 
    labs(x = "# runs", y = "estimate", title = str_c(nb_sim, " simulations")) + 
    facet_wrap(~ n_ftr)
