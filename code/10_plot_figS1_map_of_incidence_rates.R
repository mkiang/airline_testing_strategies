library(tidyverse)
library(here)
library(geofacet)
source(here("code", "mk_nytimes.R"))

state_covid <-
    read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv") %>%
    rename(fipschar = fips) %>%
    mutate(fips = as.numeric(fipschar))
pop_df <- read_csv(here("data", "nst-est2019-alldata.csv")) %>%
    select(fipschar = STATE,
           state = NAME,
           pop2019 = POPESTIMATE2019)

state_covid <- state_covid %>%
    arrange(fips, date) %>%
    left_join(pop_df) %>%
    group_by(fipschar) %>%
    mutate(new_cases = cases - lag(cases, default = 0)) %>%
    mutate(new_cases_per_100k = new_cases / pop2019 * 100000) %>%
    mutate(
        new_cases_per_100k = case_when(
            state == "Georgia" & date == as.Date("2020-10-05") ~ 0,
            new_cases_per_100k < 0 ~ 0,
            TRUE ~ new_cases_per_100k
        )
    ) %>% 
    mutate(avg_new_cases = zoo::rollmean(new_cases_per_100k, 7, fill = 0, align = "right"))

p1 <- ggplot(state_covid,
             aes(x = date, y = new_cases_per_100k)) + 
    geom_hline(
        yintercept = c(50),
        alpha = .5,
        color = "black"
    ) +
    geom_col(alpha = .5) +
    geom_line(aes(y = avg_new_cases),
              color = "red") +
    facet_geo( ~ state) +
    scale_x_date(NULL,
                 expand = c(0, 0),
                 breaks = c(as.Date("2020-04-01"),
                            as.Date("2020-07-01"),
                            as.Date("2020-10-01"),
                            as.Date("2021-01-01")),
                 labels = c("4/20", "7/20", "10/20", "1/21")
                 ) +
    scale_y_continuous("New daily reported cases per 100,000 (truncated)",
                       limits = c(0, 200),
                       expand = c(0, 0)) +
    scale_alpha_manual(NULL, values = c(.5, 1)) + 
    mk_nytimes(panel.border = element_rect(color = "grey50"))

ggsave(
    here("plots", "figS1_map_of_rolling_avg.pdf"), 
    p1,
    width = 16,
    height = 10,
    device = cairo_pdf,
    scale = .9
)
ggsave(
    here("plots", "figS1_map_of_rolling_avg.jpg"), 
    p1,
    width = 16,
    height = 10,
    dpi = 300, 
    scale = .9
)

write_csv(state_covid, here("output", "figS1_data.csv"))

## Histogram
# state_sum <- state_covid %>% 
#     group_by(state, fipschar) %>% 
#     mutate(above = (avg_new_cases >= 50) + 0) %>% 
#     summarize(days = sum(above, na.rm = TRUE), 
#               prop = mean(above, na.rm = TRUE)) %>% 
#     filter(fipschar <= 56) 
# p2 <- ggplot(state_sum, aes(x = days)) + 
#     # geom_density() + 
#     geom_histogram(binwidth = 1,
#                    color = "white") + 
#     geom_hline(yintercept = seq(0, 15, 1), 
#                color = "white",
#                alpha = 1) + 
#     scale_x_continuous("Days where 7-day average of reported new cases > 50 per 100,000",
#                        expand = c(0, 1)) + 
#     scale_y_continuous("Number of states", limits = c(0, NA), expand = c(.02, 0)) + 
#     mk_nytimes()
# 
# ggsave(
#     here("plots", "figS2_hist_rolling_avg_above_threshold.pdf"), 
#     p2,
#     width = 5,
#     height = 3.5,
#     device = cairo_pdf,
#     scale = 1.25
# )
# ggsave(
#     here("plots", "figS2_hist_rolling_avg_above_threshold.jpg"), 
#     p2,
#     width = 5,
#     height = 3.5,
#     dpi = 300, 
#     scale = 1.25
# )
# 
# sum(state_sum$prop > 0)
# sum(state_sum$days >= 14)
# median(state_sum$days)
