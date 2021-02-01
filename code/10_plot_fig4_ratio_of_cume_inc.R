## Imports ----
library(tidyverse)
library(here)
library(fs)
source(here("code", "utils.R"))
source(here("code", "mk_nytimes.R"))

plot_df <- readRDS(here("data", "summarized_results.RDS")) %>%
    filter(
        if_threshold == 0,
        time_step == 84, 
        risk_multiplier == 2,
        testing_type != "perfect_testing",
        rapid_test_multiplier == .9 |
            is.na(rapid_test_multiplier),
        sens_type != "median" | is.na(sens_type),
        prop_subclin == .3,
        symptom_adherence %in% c(0, .8),
        quarantine_adherence %in% c(0, .8),
        testing_type %in% c(
            "no_testing",
            "pcr_three_days_before",
            "pcr_three_days_before_5_day_quarantine_pcr",
            "rapid_test_same_day",
            "rapid_same_day_5_day_quarantine_pcr",
            "pcr_five_days_after"
        )
    ) %>%
    shift_time_steps() %>%
    categorize_metric() %>%
    filter(metric %in% c("cume_n_infection_daily"))

## Subset to the main analyses ----
plot_df <- bind_rows(
    ## No testing no symptom screening
    plot_df %>%
        filter(testing_type == "no_testing" &
                   symptom_adherence == 0),
    ## No quarantine 80% symptom screening adherence
    plot_df %>%
        filter(
            testing_type %in% c(
                "pcr_three_days_before",
                "rapid_test_same_day",
                "pcr_five_days_after"
            ) & 
                symptom_adherence == .8
        ),
    ## With 80% quarantine and 80% symptom screening adherence
    plot_df %>% 
        filter(
            testing_type %in% c(
                "pcr_three_days_before_5_day_quarantine_pcr",
                "rapid_same_day_5_day_quarantine_pcr"
            ) & 
                symptom_adherence == .8 & 
                quarantine_adherence == .8
        )
) 

## Add a line break to RT + PCR
levels(plot_df$testing_cat)[
    which(levels(plot_df$testing_cat) == "PCR 3 days before + 5-day quarantine")
] <-  "PCR 3 days before +\n5-day quarantine"
levels(plot_df$testing_cat)[
    which(levels(plot_df$testing_cat) == "Same-day Rapid Test + 5-day quarantine")
] <-  "Same-day Rapid Test +\n5-day quarantine"

param_grid <- expand.grid(
    destination = unique(plot_df$prob_inf),
    testing_type = unique(plot_df$testing_type),
    stringsAsFactors = FALSE
) 

holder <- vector("list", NROW(param_grid)) 

for (i in 1:NROW(holder)) {
    d_inf <- param_grid$destination[i]
    t_type <- param_grid$testing_type[i]
    
    d_cume_inf <- plot_df %>% 
        filter(testing_type == "no_testing", prob_inf == d_inf) %>%
        select(destination_prob_inf = prob_inf,
               destination_prob_inf_cat = prob_inf_cat,
               destination_inf = mean)
    
    o_cume_inf <- plot_df %>%
        filter(testing_type == t_type) %>%
        select(testing_type, 
               testing_cat, 
               origin_prob_inf = prob_inf, 
               origin_prob_inf_cat = prob_inf_cat,
               origin_inf = mean)
    
    holder[[i]] <- o_cume_inf %>% 
        bind_cols(d_cume_inf) %>% 
        mutate(ratio = origin_inf / destination_inf)
}

holder <- bind_rows(holder)

holder <- holder %>%
    mutate(ratio_trunc = case_when(ratio <= 1 ~ 1,
                                   ratio >= 15 ~ 15,
                                   TRUE ~ ratio))

levels(holder$origin_prob_inf_cat) <- gsub(" per 100,000", "", levels(holder$origin_prob_inf_cat))
levels(holder$destination_prob_inf_cat) <- gsub(" per 100,000", "", levels(holder$destination_prob_inf_cat))

p1 <- ggplot(
    holder %>%
        filter(
            origin_prob_inf >= destination_prob_inf,
            between(origin_prob_inf, 20 / 100000, 500 / 100000),
            between(destination_prob_inf, 20 / 100000, 500 / 100000)
        ),
    aes(
        x = as.factor(destination_prob_inf_cat),
        y = as.factor(origin_prob_inf_cat),
        fill = ratio_trunc
    )
) +
    geom_tile(color = "white") +
    facet_wrap( ~ testing_cat) +
    scale_y_discrete("Daily probability of infection at origin per 100,000",
                     expand = c(0, 0)) +
    scale_x_discrete("Daily probability of infection at destination per 100,000",
                     expand = c(0, 0)) +
    scale_fill_distiller(
        "Ratio of cumulative\ninfectious days",
        palette = "YlOrRd",
        direction = 1,
        trans = "log2",
        breaks = c(1, 2, 4, 8, 15),
        labels = c(1, 2, 4, 8, "15+"),
        guide = guide_colorbar(barwidth = unit(.5, "cm"),
                               barheight = unit(6, "cm"))
    ) +
    mk_nytimes(
        legend.position = "right",
        # axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank()
    ) +
    coord_equal()

ggsave(
    "./plots/fig4_ratio_of_cume_inc.pdf",
    p1,
    device = cairo_pdf,
    width = 8,
    height = 5,
    scale = 1.1
)

ggsave(
    "./plots/fig4_ratio_of_cume_inc.jpg",
    p1,
    dpi = 300, 
    width = 8,
    height = 5,
    scale = 1.1
)

write_csv(
    holder %>%
        filter(
            origin_prob_inf >= destination_prob_inf,
            between(origin_prob_inf, 20 / 100000, 500 / 100000),
            between(destination_prob_inf, 20 / 100000, 500 / 100000)
        ),
    "./output/fig4_data.csv"
)

# p1 <- ggplot(
#     holder %>%
#         filter(
#             origin_prob_inf >= destination_prob_inf,
#             between(origin_prob_inf, 10 / 100000, 500 / 100000),
#             between(destination_prob_inf, 10 / 100000, 500 / 100000)
#         ),
#     aes(
#         x = as.factor(destination_prob_inf_cat),
#         y = as.factor(origin_prob_inf_cat),
#         fill = ratio_trunc
#     )
# ) +
#     geom_tile(color = "white") +
#     facet_wrap( ~ testing_cat) +
#     scale_y_discrete("Daily probability of infection at origin per 100,000",
#                      expand = c(0, 0)) +
#     scale_x_discrete("Daily probability of infection at destination per 100,000",
#                      expand = c(0, 0)) +
#     scale_fill_distiller(
#         "Ratio of cumulative\ninfectious days",
#         palette = "YlOrRd",
#         direction = 1,
#         trans = "log2",
#         breaks = c(1, 2, 4, 8, 15),
#         labels = c(1, 2, 4, 8, "15+"),
#         guide = guide_colorbar(barwidth = unit(.5, "cm"),
#                                barheight = unit(6, "cm"))
#     ) +
#     mk_nytimes(
#         legend.position = "right",
#         # axis.text.x = element_text(angle = 90, hjust = 1),
#         panel.grid.major = element_blank()
#     ) +
#     coord_equal()
# 
# ggsave(
#     "./plots/fig4a_ratio_of_cume_inc.pdf",
#     p1,
#     device = cairo_pdf,
#     width = 8,
#     height = 5,
#     scale = 1.1
# )
# 
# ggsave(
#     "./plots/fig4a_ratio_of_cume_inc.jpg",
#     p1,
#     dpi = 300, 
#     width = 8,
#     height = 5,
#     scale = 1.1
# )
# 
# p1 <- ggplot(
#     holder %>%
#         filter(
#             origin_prob_inf >= destination_prob_inf,
#             between(origin_prob_inf, 20 / 100000, 250 / 100000),
#             between(destination_prob_inf, 20 / 100000, 250 / 100000)
#         ),
#     aes(
#         x = as.factor(destination_prob_inf_cat),
#         y = as.factor(origin_prob_inf_cat),
#         fill = ratio_trunc
#     )
# ) +
#     geom_tile(color = "white") +
#     facet_wrap( ~ testing_cat) +
#     scale_y_discrete("Daily probability of infection at origin per 100,000",
#                      expand = c(0, 0)) +
#     scale_x_discrete("Daily probability of infection at destination per 100,000",
#                      expand = c(0, 0)) +
#     scale_fill_distiller(
#         "Ratio of cumulative\ninfectious days",
#         palette = "YlOrRd",
#         direction = 1,
#         trans = "log2",
#         breaks = c(1, 2, 4, 8, 15),
#         labels = c(1, 2, 4, 8, "15+"),
#         guide = guide_colorbar(barwidth = unit(.5, "cm"),
#                                barheight = unit(6, "cm"))
#     ) +
#     mk_nytimes(
#         legend.position = "right",
#         # axis.text.x = element_text(angle = 90, hjust = 1),
#         panel.grid.major = element_blank()
#     ) +
#     coord_equal()
# 
# ggsave(
#     "./plots/fig4c_ratio_of_cume_inc.pdf",
#     p1,
#     device = cairo_pdf,
#     width = 8,
#     height = 5,
#     scale = 1.1
# )
# 
# ggsave(
#     "./plots/fig4c_ratio_of_cume_inc.jpg",
#     p1,
#     dpi = 300, 
#     width = 8,
#     height = 5,
#     scale = 1.1
# )
