## Imports ----
library(tidyverse)
library(here)
library(fs)
source(here::here("code", "utils.R"))
source(here::here("code", "mk_nytimes.R"))

## Constants ----
PROB_INF <- config::get()$primary_prob_inf / 1000000

plot_df <-
    readRDS(here::here("data", "summarized_testing_results.RDS"))  %>%
    dplyr::filter(metric %in% c("n_active_infected_day_of_flight")) %>%
    dplyr::filter(
        if_threshold == 0,
        prob_inf == PROB_INF,
        risk_multiplier == 2,
        testing_type != "perfect_testing",
        rapid_test_multiplier == .9 |
            is.na(rapid_test_multiplier),
        sens_type != "median" | is.na(sens_type),
        prop_subclin == .3,
        testing_type %in% c(
            "no_testing",
            "pcr_three_days_before",
            # "pcr_three_days_before_5_day_quarantine_pcr",
            "rapid_test_same_day",
            # "rapid_same_day_5_day_quarantine_pcr",
            "pcr_five_days_after"
        )
    ) %>%
    shift_time_steps() %>%
    categorize_metric()

plot_df <- plot_df %>%
    dplyr::mutate(symptom_adherence_cat = factor(
        symptom_adherence,
        levels = seq(0, 1, .2),
        labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
        ordered = TRUE
    ))

p1 <- ggplot2::ggplot(plot_df,
                      ggplot2::aes(
                          x = symptom_adherence_cat,
                          y = mean,
                          ymax = p975,
                          ymin = p025
                      )) +
    ggplot2::geom_point(alpha = 1) +
    ggplot2::geom_errorbar(width = .1) +
    ggplot2::facet_wrap( ~ testing_cat, nrow = 1) +
    mk_nytimes(panel.border = ggplot2::element_rect(color = "grey70")) +
    ggplot2::scale_x_discrete("Percent of passengers adherening to symptomatic self-isolation") +
    ggplot2::scale_y_continuous("Number of actively infectious passengers on day of flight")

ggplot2::ggsave(
    "./plots/figS9_active_inf_dof.pdf",
    p1,
    device = grDevices::cairo_pdf,
    width = 7,
    height = 4,
    scale = 1.2
)

ggplot2::ggsave(
    "./plots/figS9_active_inf_dof.jpg",
    p1,
    dpi = 300,
    width = 7,
    height = 4,
    scale = 1.2
)

readr::write_csv(plot_df,
                 "./output/figS9_data.csv")
