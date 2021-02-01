## Imports ----
library(tidyverse)
library(here)
library(fs)
library(patchwork)
source(here::here("code", "utils.R"))
source(here::here("code", "mk_nytimes.R"))

## Constants ----
plot_df <-
    readRDS(here::here("data", "summarized_testing_results.RDS")) %>%
    dplyr::filter(
        if_threshold == 0,
        risk_multiplier == 2,
        prop_subclin == .3,
        sens_type == "upper" | is.na(sens_type),
        symptom_adherence %in% c(0, .8),
        testing_type %in% c("no_testing", "pcr_three_days_before", "rapid_test_same_day"),
        metric == "n_active_infected_day_of_flight",
        rapid_test_multiplier == .9 |
            is.na(rapid_test_multiplier)
    ) %>%
    categorize_metric()

plot_df <- dplyr::bind_rows(
    plot_df %>%
        dplyr::filter(testing_type == "no_testing" &
                          symptom_adherence == 0),
    plot_df %>%
        dplyr::filter(testing_type != "no_testing" &
                          symptom_adherence == .8)
)

p1 <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
        x = prob_inf,
        y = mean,
        ymin = p025,
        ymax = p975,
        color = testing_cat,
        fill = testing_cat,
        group = testing_cat
    )
) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = .2,
                         color = NA) +
    ggplot2::scale_x_continuous(
        "Daily probability of infection per 100,000",
        expand = c(0, 0),
        breaks = unique(plot_df$prob_inf)[-2],
        labels = function(x)
            as.numeric(x) * 100000
    ) +
    ggplot2::scale_y_continuous("Mean (95% UI)",
                                expand = c(0, 0),
                                limits = c(0, NA)) +
    ggplot2::scale_color_manual("Testing strategy",
                                values = c("black",
                                           RColorBrewer::brewer.pal(name = "Dark2", n = 5)[1:2])) +
    ggplot2::scale_fill_manual("Testing strategy",
                               values = c("black",
                                          RColorBrewer::brewer.pal(name = "Dark2", n = 5)[1:2])) +
    mk_nytimes(
        panel.border = ggplot2::element_rect(color = "grey30"),
        legend.background = ggplot2::element_rect(color = NA,
                                                  fill = ggplot2::alpha("white", .75)),
        legend.position = c(.01, .99),
        legend.justification = c(0, 1),
        axis.text.x = ggplot2::element_text(hjust = c(0, .5, .5, .5, .5, .5, 1))
    )

ggplot2::ggsave(
    "./plots/fig3_active_inf_day_of_travel.pdf",
    p1,
    device = grDevices::cairo_pdf,
    width = 5,
    height = 3,
    scale = 1.2
)

ggplot2::ggsave(
    "./plots/fig3_active_inf_day_of_travel.jpg",
    p1,
    dpi = 300,
    width = 5,
    height = 3,
    scale = 1.2
)

readr::write_csv(plot_df,
                 "./output/fig3_data.csv")
