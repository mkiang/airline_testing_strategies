library(tidyverse)
library(here)
library(fs)
library(foreach)
library(doParallel)
library(patchwork)
source(here("code", "mk_nytimes.R"))

## Read in all simulated populations and extract the random draws ----
if (!file_exists(here("data", "observed_parameter_draws.RDS"))) {
    sim_files <- dir_ls(here("intermediate_files", "simulation_states"), recurse = TRUE, glob = "*.RDS")
    
    doParallel::registerDoParallel(cores = 15)
    sim_parameters <- foreach::foreach(i = 1:NROW(sim_files)) %dopar% {
        x <- readRDS(sim_files[i])
        tibble(
            p_spec = x$p_spec,
            n_incubation_days = x$days_incubation,
            n_presxinf_days = as.numeric(names(x$days_incubation)),
            n_sxinf_days = as.numeric(names(x$days_symptomatic)),
            n_symptomatic_days = x$days_symptomatic
        ) %>% 
            mutate(n_inf_days = n_sxinf_days + n_presxinf_days,
                   n_presxnoinf_days = n_incubation_days - n_presxinf_days,
                   n_sxpostinf_days = n_symptomatic_days - n_sxinf_days)
    }
    sim_parameters <- bind_rows(sim_parameters)
    
    saveRDS(sim_parameters, here("data", "observed_parameter_draws.RDS"))
} else {
    sim_parameters <- readRDS(here("data", "observed_parameter_draws.RDS"))
}

p1 <- ggplot(sim_parameters,
                   aes(x = p_spec)) +
    geom_density(fill = "black", 
                 alpha = .2, 
                 color = "black") + 
    scale_x_continuous("Specificity", 
                       expand = c(.025, 0)) + 
    scale_y_continuous("Density", expand = c(.025, 0)) + 
    mk_nytimes()

p2 <- ggplot(sim_parameters,
       aes(x = n_presxnoinf_days)) +
    geom_histogram(binwidth = 1, 
                   fill = "black", 
                   alpha = .2, 
                   color = "black") + 
    scale_x_continuous("Incubation days", 
                       expand = c(.025, 0)) + 
    scale_y_continuous("Frequency", expand = c(.025, 0)) + 
    mk_nytimes()

p3 <- ggplot(sim_parameters,
       aes(x = n_presxinf_days)) +
    geom_histogram(binwidth = 1, 
                   fill = "black", 
                   alpha = .2, 
                   color = "black") + 
    scale_x_continuous("Pre-symptomatic infectious days", 
                       expand = c(.025, 0)) + 
    scale_y_continuous(NULL, expand = c(.025, 0)) + 
    mk_nytimes()

p4 <- ggplot(sim_parameters,
       aes(x = n_sxinf_days)) +
    geom_histogram(binwidth = 1, 
                   fill = "black", 
                   alpha = .2, 
                   color = "black") + 
    scale_x_continuous("Symptomatic infectious days", 
                       expand = c(.025, 0)) + 
    scale_y_continuous(NULL, expand = c(.025, 0)) + 
    mk_nytimes()

p5 <- ggplot(sim_parameters,
       aes(x = n_inf_days)) +
    geom_histogram(binwidth = 1, 
                   fill = "black", 
                   alpha = .2, 
                   color = "black") + 
    scale_x_continuous("Total infectious days", 
                       expand = c(.025, 0)) + 
    scale_y_continuous(NULL, expand = c(.025, 0)) + 
    mk_nytimes()
p_all <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 1)

ggsave(here("plots", "figS3_distribution_of_simulation_draws.pdf"), 
       p_all,
       width = 12.5,
       height = 3,
       scale = 1,
       device = cairo_pdf)
ggsave(here("plots", "figS3_distribution_of_simulation_draws.jpg"), 
       p_all,
       width = 12.5,
       height = 3,
       scale = 1,
       dpi = 300)
