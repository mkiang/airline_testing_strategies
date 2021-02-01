library(tidyverse)
library(here)
library(geofacet)
source(here::here("code", "mk_nytimes.R"))

state_covid <-
    readr::read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv") %>%
    dplyr::rename(fipschar = fips) %>%
    dplyr::mutate(fips = as.numeric(fipschar))
pop_df <-
    readr::read_csv(here::here("data", "nst-est2019-alldata.csv")) %>%
    dplyr::select(fipschar = STATE,
                  state = NAME,
                  pop2019 = POPESTIMATE2019)

state_covid <- state_covid %>%
    dplyr::arrange(fips, date) %>%
    dplyr::left_join(pop_df) %>%
    dplyr::group_by(fipschar) %>%
    dplyr::mutate(new_cases = cases - dplyr::lag(cases, default = 0)) %>%
    dplyr::mutate(new_cases_per_100k = new_cases / pop2019 * 100000) %>%
    dplyr::mutate(
        new_cases_per_100k = dplyr::case_when(
            state == "Georgia" & date == as.Date("2020-10-05") ~ 0,
            new_cases_per_100k < 0 ~ 0,
            TRUE ~ new_cases_per_100k
        )
    ) %>%
    dplyr::mutate(avg_new_cases = zoo::rollmean(new_cases_per_100k, 7, fill = 0, align = "right"))

p1 <- ggplot2::ggplot(state_covid,
                      ggplot2::aes(x = date, y = new_cases_per_100k)) +
    ggplot2::geom_hline(yintercept = c(50),
                        alpha = .5,
                        color = "black") +
    ggplot2::geom_col(alpha = .5) +
    ggplot2::geom_line(ggplot2::aes(y = avg_new_cases),
                       color = "red") +
    geofacet::facet_geo(~ state) +
    ggplot2::scale_x_date(
        NULL,
        expand = c(0, 0),
        breaks = c(
            as.Date("2020-04-01"),
            as.Date("2020-07-01"),
            as.Date("2020-10-01"),
            as.Date("2021-01-01")
        ),
        labels = c("Apr", "Jul", "Oct", "Jan")
    ) +
    ggplot2::scale_y_continuous(
        "New daily reported cases per 100,000 (truncated)",
        limits = c(0, 200),
        expand = c(0, 0)
    ) +
    ggplot2::scale_alpha_manual(NULL, values = c(.5, 1)) +
    mk_nytimes(panel.border = ggplot2::element_rect(color = "grey50"))

ggplot2::ggsave(
    here::here("plots", "figS1_map_of_rolling_avg.pdf"),
    p1,
    width = 16,
    height = 10,
    device = grDevices::cairo_pdf,
    scale = .9
)
ggplot2::ggsave(
    here::here("plots", "figS1_map_of_rolling_avg.jpg"),
    p1,
    width = 16,
    height = 10,
    dpi = 300,
    scale = .9
)

readr::write_csv(state_covid, here::here("output", "figS1_data.csv"))
