## Imports ----
library(tidyverse)
library(here)
library(fs)
source(here("code", "utils.R"))
source(here("code", "mk_nytimes.R"))

plot_df <- readRDS(here("data", "summarized_testing_results.RDS"))

plot_df <- plot_df %>%
    filter(
        risk_multiplier == 2,
        testing_type != "perfect_testing", 
        metric == "ratio_false_true",
        testing_type != "no_testing", 
        rapid_test_multiplier == .9 | is.na(rapid_test_multiplier),
        prob_inf > 5 / 1000000,
        prop_subclin == .3, 
        sens_type == "upper", 
        if_threshold == 0,
        symptom_adherence == .8
    ) %>%
    categorize_metric()

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
        x = as.factor(prob_inf),
        y = mean,
        ymin = p025,
        ymax = p975,
        group = testing_cat,
        color = testing_cat
    )
) +
    geom_hline(yintercept = 1, alpha = .5) +
    geom_point(position = position_dodge(width = .5)) +
    geom_errorbar(width = .1, position = position_dodge(width = .5)) +
    scale_x_discrete(
        "Daily probability of infection per 100,000",
        expand = c(0, 0),
        labels = function(x)
            formatC(as.numeric(x) * 100000, big.mark = ",")
    ) +
    scale_y_continuous(
        "Mean (95% UI) ratio of False/True positive results",
        expand = c(0, 0),
        breaks = c(
            0,
            # 1 / 200,
            1 / 100,
            1 / 50,
            1 / 25,
            1 / 10,
            1 / 5,
            1 / 2,
            1,
            2,
            5,
            10,
            25,
            50,
            100,
            200
        ),
        trans = "log"
    ) +
    scale_color_brewer("Testing strategy",
                       palette = "Dark2",
                       drop = TRUE) +
    mk_nytimes(
        panel.border = element_rect(color = "grey30"),
        legend.background = element_rect(color = NA,
                                         fill = alpha("white", .75)),
        legend.position = "right" # ,
        # legend.justification = c(1, 1)
    )

ggsave(
    "./plots/fig2_falsetrue_ratio.pdf",
    p1,
    device = cairo_pdf,
    width = 9,
    height = 4.5,
    scale = 1.1
)

ggsave(
    "./plots/fig2_falsetrue_ratio.jpg",
    p1,
    dpi = 300, 
    width = 9,
    height = 4,
    scale = 1.1
)

write_csv(
    plot_df,
    "./output/fig2_data.csv"
)
