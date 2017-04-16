# Load required packages --------------------------------------------------
library(tidyverse)
library(stringr)
library(broom)


# Simulation data ---------------------------------------------------------
nb_sim <- 100
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
    gather(estimate, ave_sim, ave_sd, ave_mad) %>% 
    mutate(estimate = ifelse(str_detect(estimate, "sd"), "SD", "MAD"))

sd_sum <- df_sum %>% 
    select(-ave_sd, -ave_mad) %>% 
    gather(estimate, sd_sim, std_sd, std_mad) %>% 
    mutate(estimate = ifelse(str_detect(estimate, "sd"), "SD", "MAD"))

df_sum2 <- left_join(ave_sum, sd_sum)

df_sum2 %>% ggplot(aes(size, ave_sim, colour = estimate)) + 
    geom_point()


