## Imports ----
library(tidyverse)
library(here)
library(fs)
source(here::here("code", "utils.R"))
source(here::here("code", "mk_nytimes.R"))

## Constants ----
PROB_INF <- config::get()$primary_prob_inf / 1000000

plot_df <- readRDS(here::here("data", "summarized_results.RDS")) %>%
    dplyr::filter(
        if_threshold == 0,
        prob_inf == PROB_INF,
        risk_multiplier == 2,
        testing_type != "perfect_testing",
        rapid_test_multiplier == .9 |
            is.na(rapid_test_multiplier),
        prop_subclin == .3,
        sens_type != "median" | is.na(sens_type),
        symptom_adherence %in% c(0),
        quarantine_adherence %in% c(0),
        testing_type %in% c(
            "pcr_three_days_before_5_day_quarantine_pcr",
            "pcr_three_days_before_7_day_quarantine_pcr",
            "pcr_three_days_before_14_day_quarantine_pcr",
            "rapid_same_day_5_day_quarantine_pcr",
            "rapid_same_day_7_day_quarantine_pcr",
            "rapid_same_day_14_day_quarantine_pcr"
        )
    ) %>%
    shift_time_steps() %>%
    categorize_metric() %>%
    dplyr::filter(metric %in% c("cume_n_infection_daily")) %>%
    categorize_testing_types()

p1 <- ggplot2::ggplot(plot_df,
                      ggplot2::aes(
                          x = relative_time,
                          y = mean,
                          ymin = p025,
                          ymax = p975,
                      )) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "dashed",
                        alpha = .8) +
    ggplot2::geom_ribbon(color = NA, alpha = .2) +
    ggplot2::geom_line(alpha = .8) +
    ggplot2::facet_wrap( ~ testing_cat, nrow = 2) +
    ggplot2::scale_x_continuous(
        "Time relative to flight, days",
        expand = c(0, 0),
        breaks = c(-3, 0, 7, 14)
    ) +
    ggplot2::scale_y_continuous("Mean (95% UI) cumulative infectious days",
                                expand = c(0, 0)) +
    ggplot2::scale_color_brewer("Percent of infections\nthat are subclinical",
                                palette = "Set1") +
    ggplot2::scale_fill_brewer("Percent of infections\nthat are subclinical",
                               palette = "Set1") +
    mk_nytimes(
        panel.border = ggplot2::element_rect(color = "grey30"),
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(hjust = c(0, .5, .5, 1))
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1),
                                                   title.position = "top"))

ggplot2::ggsave(
    "./plots/figS13_by_quarantine_duration.pdf",
    p1,
    device = grDevices::cairo_pdf,
    width = 8,
    height = 5,
    scale = 1.1
)

ggplot2::ggsave(
    "./plots/figS13_by_quarantine_duration.jpg",
    p1,
    dpi = 300,
    width = 8,
    height = 5,
    scale = 1.1
)

readr::write_csv(plot_df,
                 "./output/figS13_data.csv")
