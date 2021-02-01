library(tidyverse)
library(here)
library(fs)
library(foreach)
library(doParallel)
library(patchwork)
source(here::here("code", "mk_nytimes.R"))

## Read in all simulated populations and extract the random draws ----
if (!fs::file_exists(here::here("data", "observed_parameter_draws.RDS"))) {
    ## Get a list of all the primary simulation state files
    sim_params <- vector("list", 50)
    for (i in 1:NROW(sim_params)) {
        sim_params[[i]] <- return_parameter_grid(new_only = FALSE) %>%
            dplyr::mutate(rep = i)
    }
    sim_params <- dplyr::bind_rows(sim_params)
    sim_files <- unique(
        return_sim_state_file(
            sim_params$round,
            sim_params$rep,
            sim_params$prob_inf,
            sim_params$prop_subclin
        )
    )
    
    doParallel::registerDoParallel(cores = 15)
    sim_parameters <-
        foreach::foreach(i = 1:NROW(sim_files)) %dopar% {
            x <- readRDS(sim_files[i])
            dplyr::tibble(
                p_spec = x$p_spec,
                n_incubation_days = x$days_incubation,
                n_presxinf_days = as.numeric(names(x$days_incubation)),
                n_sxinf_days = as.numeric(names(x$days_symptomatic)),
                n_symptomatic_days = x$days_symptomatic
            ) %>%
                dplyr::mutate(
                    n_inf_days = n_sxinf_days + n_presxinf_days,
                    n_presxnoinf_days = n_incubation_days - n_presxinf_days,
                    n_sxpostinf_days = n_symptomatic_days - n_sxinf_days
                )
        }
    sim_parameters <- dplyr::bind_rows(sim_parameters)
    
    saveRDS(sim_parameters,
            here::here("data", "observed_parameter_draws.RDS"))
} else {
    sim_parameters <-
        readRDS(here::here("data", "observed_parameter_draws.RDS"))
}

p1 <- ggplot2::ggplot(sim_parameters,
                      ggplot2::aes(x = p_spec)) +
    ggplot2::geom_density(fill = "black",
                          alpha = .2,
                          color = "black") +
    ggplot2::scale_x_continuous("Specificity",
                                expand = c(.025, 0)) +
    ggplot2::scale_y_continuous("Density", expand = c(.025, 0)) +
    mk_nytimes()

p2 <- ggplot2::ggplot(sim_parameters,
                      ggplot2::aes(x = n_presxnoinf_days)) +
    ggplot2::geom_histogram(
        binwidth = 1,
        fill = "black",
        alpha = .2,
        color = "black"
    ) +
    ggplot2::scale_x_continuous("Incubation days",
                                expand = c(.025, 0)) +
    ggplot2::scale_y_continuous("Frequency", expand = c(.025, 0)) +
    mk_nytimes()

p3 <- ggplot2::ggplot(sim_parameters,
                      ggplot2::aes(x = n_presxinf_days)) +
    ggplot2::geom_histogram(
        binwidth = 1,
        fill = "black",
        alpha = .2,
        color = "black"
    ) +
    ggplot2::scale_x_continuous("Pre-symptomatic infectious days",
                                expand = c(.025, 0)) +
    ggplot2::scale_y_continuous(NULL, expand = c(.025, 0)) +
    mk_nytimes()

p4 <- ggplot2::ggplot(sim_parameters,
                      ggplot2::aes(x = n_sxinf_days)) +
    ggplot2::geom_histogram(
        binwidth = 1,
        fill = "black",
        alpha = .2,
        color = "black"
    ) +
    ggplot2::scale_x_continuous("Symptomatic infectious days",
                                expand = c(.025, 0)) +
    ggplot2::scale_y_continuous(NULL, expand = c(.025, 0)) +
    mk_nytimes()

p5 <- ggplot2::ggplot(sim_parameters,
                      ggplot2::aes(x = n_inf_days)) +
    ggplot2::geom_histogram(
        binwidth = 1,
        fill = "black",
        alpha = .2,
        color = "black"
    ) +
    ggplot2::scale_x_continuous("Total infectious days",
                                expand = c(.025, 0)) +
    ggplot2::scale_y_continuous(NULL, expand = c(.025, 0)) +
    mk_nytimes()
p_all <- p1 + p2 + p3 + p4 + p5 + patchwork::plot_layout(nrow = 1)

ggplot2::ggsave(
    here::here("plots", "figS2_distribution_of_simulation_draws.pdf"),
    p_all,
    width = 12.5,
    height = 3,
    scale = 1,
    device = grDevices::cairo_pdf
)
ggplot2::ggsave(
    here::here("plots", "figS2_distribution_of_simulation_draws.jpg"),
    p_all,
    width = 12.5,
    height = 3,
    scale = 1,
    dpi = 300
)
