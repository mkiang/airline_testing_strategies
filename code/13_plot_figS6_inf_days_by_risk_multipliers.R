## Imports ----
library(tidyverse)
library(here)
library(fs)
source(here("code", "utils.R"))
source(here("code", "mk_nytimes.R"))

## Constants ----
PROB_INF <- config::get()$primary_prob_inf / 1000000

plot_df <- readRDS(here("data", "summarized_results.RDS")) %>% 
    filter(if_threshold == 0, 
           testing_type != "perfect_testing", 
           prob_inf == PROB_INF,  
           rapid_test_multiplier == .9 | is.na(rapid_test_multiplier)) %>% 
    shift_time_steps() %>% 
    categorize_metric()

## Add a line break to RT + PCR
levels(plot_df$testing_cat)[7] <- "PCR 3 days before +\nSame-day RT" 
levels(plot_df$testing_cat)[8] <- "Same-day RT + 5-day\nquarantine + PCR"
levels(plot_df$testing_cat)[9] <- "PCR 3 days before +\n5-day quarantine + PCR"

p1 <- ggplot(
    plot_df %>%
        filter(
            metric %in% c("cume_n_infection_daily")
        ),
    aes(
        x = relative_time,
        y = mean,
        ymin = p025,
        ymax = p975, 
        group = symptom_cat,
        color = symptom_cat,
        fill = symptom_cat
    )
) +
    geom_vline(xintercept = 0,
               linetype = "dashed",
               alpha = .8) +
    geom_ribbon(color = NA, alpha = .2) +  
    geom_line(alpha = .8) +
    facet_grid(risk_multi_cat ~ testing_cat,
               scales = "free") +
    scale_x_continuous(
        "Time relative to flight, days",
        expand = c(0, 0),
        breaks = c(-3, 0, 7, 14)
    ) +
    scale_y_continuous("Mean (95% UI) cumulative infectious days",
                       expand = c(0, 0)) +
    scale_color_brewer(NULL,
                       palette = "Set1") +
    scale_fill_brewer(NULL,
                      palette = "Set1") +
    mk_nytimes(
        panel.border = element_rect(color = "grey30"),
        legend.position = "bottom",
        axis.text.x = element_text(hjust = c(0, .5, .5, 1))
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),
                                 title.position = "top")) 

ggsave(
    "./plots/figS6_cume_infdays_by_risk_multi.pdf",
    p1,
    device = cairo_pdf,
    width = 9,
    height = 9,
    scale = 1.1
)

ggsave(
    "./plots/figS6_cume_infdays_by_risk_multi.jpg",
    p1,
    dpi = 300, 
    width = 9,
    height = 9,
    scale = 1.1
)

write_csv(
    plot_df %>%
        filter(
            metric %in% c("cume_n_infection_daily")
        ),
    "./output/figS6_data.csv"
)
