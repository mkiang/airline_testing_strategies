## Imports ----
library(tidyverse)
library(here)
library(fs)
library(patchwork)
source(here("code", "utils.R"))
source(here("code", "mk_nytimes.R"))

## Constants ----
plot_df <- readRDS(here("data", "summarized_testing_results.RDS")) %>%
    filter(
        if_threshold == 0,
        risk_multiplier == 2,
        testing_type %in% c("no_testing", "pcr_three_days_before", "rapid_test_same_day"),
        rapid_test_multiplier == .9 |
            is.na(rapid_test_multiplier)
    ) %>% 
    categorize_metric()

plot_df <- bind_rows(
    plot_df %>% 
        filter(testing_type == "no_testing" & !symptom_screening), 
    plot_df %>% 
        filter(testing_type != "no_testing" & symptom_screening)
)

## Add a line break to RT + PCR
levels(plot_df$testing_cat)[7] <- "PCR 3 days before +\nSame-day RT" 
levels(plot_df$testing_cat)[8] <- "Same-day RT + 5-day\nquarantine + PCR"
levels(plot_df$testing_cat)[9] <- "PCR 3 days before +\n5-day quarantine + PCR" 

p1 <- ggplot(
    plot_df %>%
        filter(metric %in% c("n_active_infected_day_of_flight")),
    aes(
        x = prob_inf,
        y = mean,
        ymin = p025,
        ymax = p975,
        color = testing_cat,
        fill = testing_cat,
        group = testing_cat
    )
) +
    geom_line() +
    geom_ribbon(alpha = .2,
                color = NA) +
    scale_x_continuous(
        "Daily probability of infection per 100,000",
        expand = c(0, 0),
        breaks = c(1, seq(50, 250, 50)) / 100000, 
        labels = function(x)
            as.numeric(x) * 100000
    ) +
    scale_y_continuous("Mean (95% UI)",
                       expand = c(0, 0),
                       limits = c(0, NA)) +
    scale_color_manual("Testing strategy",
                       values = c("black",
                                  RColorBrewer::brewer.pal(name = "Dark2", n = 5)[1:2])) +
    scale_fill_manual("Testing strategy",
                      values = c("black",
                                 RColorBrewer::brewer.pal(name = "Dark2", n = 5)[1:2])) +
    mk_nytimes(
        panel.border = element_rect(color = "grey30"),
        legend.background = element_rect(color = NA,
                                         fill = alpha("white", .75)),
        legend.position = c(.01, .99),
        legend.justification = c(0, 1),
        axis.text.x = element_text(hjust = c(0, .5, .5, .5, .5, 1))
    )

ggsave(
    "./plots/fig3_active_inf_day_of_travel.pdf",
    p1,
    device = cairo_pdf,
    width = 5,
    height = 3,
    scale = 1.2
)

ggsave(
    "./plots/fig3_active_inf_day_of_travel.jpg",
    p1,
    dpi = 300, 
    width = 5,
    height = 3,
    scale = 1.2
)

write_csv(
    plot_df %>%
        filter(metric %in% c("n_active_infected_day_of_flight"), 
               symptom_screening),
    "./output/fig3_data.csv"
)
