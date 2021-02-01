## Imports ----
library(tidyverse)
library(here)
library(fs)
source(here("code", "utils.R"))
source(here("code", "mk_nytimes.R"))

## Constants ----
PROB_INF <- config::get()$primary_prob_inf / 1000000

plot_df <- readRDS(here("data", "summarized_results.RDS")) %>%
    filter(
        if_threshold == 0,
        prob_inf == PROB_INF,
        risk_multiplier == 2,
        testing_type != "perfect_testing",
        rapid_test_multiplier == .9 |
            is.na(rapid_test_multiplier),
        sens_type != "median" | is.na(sens_type),
        prop_subclin == .3,
        symptom_adherence ==.8,
        testing_type %in% c(
            "pcr_three_days_before_5_day_quarantine_pcr",
            "rapid_same_day_5_day_quarantine_pcr"
        )
    ) %>%
    shift_time_steps() %>%
    categorize_metric() %>%
    filter(metric %in% c("cume_n_infection_daily")) %>% 
    mutate(quarantine_adherence_cat = factor(
        quarantine_adherence,
        levels = seq(0, 1, .2),
        labels = c("0%", "20%", "40%", "60%", "80% (baseline)", "100%"),
        ordered = TRUE
    ))

## Add a line break to RT + PCR
levels(plot_df$testing_cat)[
    which(levels(plot_df$testing_cat) == "PCR 3 days before + 5-day quarantine")
] <-  "PCR 3 days before +\n5-day quarantine"
levels(plot_df$testing_cat)[
    which(levels(plot_df$testing_cat) == "Same-day Rapid Test + 5-day quarantine")
] <-  "Same-day Rapid Test +\n5-day quarantine"

p1 <- ggplot(
    plot_df,
    aes(
        x = relative_time,
        y = mean,
        ymin = p025,
        ymax = p975
    )
) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               alpha = .8) +
    geom_ribbon(color = NA, alpha = .2) +  
    geom_line(alpha = .8) +
    facet_grid(quarantine_adherence_cat ~ testing_cat) +
    scale_x_continuous(
        "Time relative to flight, days",
        expand = c(0, 0),
        breaks = c(-3, 0, 7, 14)
    ) +
    scale_y_continuous("Mean (95% UI) cumulative infectious days",
                       expand = c(0, 0)) +
    scale_color_brewer("Day-specific test sensitivity",
                       labels = c("Median", "Upper bound"), 
                       palette = "Set1") +
    scale_fill_brewer("Day-specific test sensitivity",
                      labels = c("Median", "Upper bound"), 
                      palette = "Set1") +
    mk_nytimes(
        panel.border = element_rect(color = "grey30"),
        legend.position = "bottom",
        axis.text.x = element_text(hjust = c(0, .5, .5, 1))
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),
                                 title.position = "top")) 

ggsave(
    "./plots/figS10_cume_inf_quarantine_adherence.pdf",
    p1,
    device = cairo_pdf,
    width = 4,
    height = 8,
    scale = 1.1
)

ggsave(
    "./plots/figS10_cume_inf_quarantine_adherence.jpg",
    p1,
    dpi = 300,
    width = 4,
    height = 8,
    scale = 1.1
)

write_csv(
    plot_df,
    "./output/figS10_data.csv"
)