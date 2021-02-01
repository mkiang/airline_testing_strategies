library(tidyverse)
library(here)
source(here::here("code", "utils.R"))
source(here::here("code", "mk_nytimes.R"))

plot_df <- sens_by_day_data() %>%
    dplyr::mutate(day = 1:dplyr::n(),
                  rt = .9 * upper)

p2 <- ggplot2::ggplot(plot_df,
                      ggplot2::aes(
                          x = day,
                          ymin = lower,
                          ymax = upper,
                          y = median
                      )) +
    ggplot2::geom_ribbon(alpha = .15) +
    ggplot2::geom_point(alpha = .4) +
    ggplot2::geom_line(
        data = plot_df,
        ggplot2::aes(x = day, y = upper),
        color = "black",
        alpha = 1,
        size = .75
    ) +
    ggplot2::scale_x_continuous("Time from exposure, days",
                                expand = c(.025, 0)) +
    ggplot2::scale_y_continuous(
        "Test sensitivity, %",
        expand = c(.025, 0),
        labels = function(x)
            round(x * 100)
    ) +
    mk_nytimes()

ggplot2::ggsave(
    here::here("plots", "figS3_test_sens.pdf"),
    p2,
    width = 5,
    height = 3.5,
    scale = 1,
    device = grDevices::cairo_pdf
)

ggplot2::ggsave(
    here::here("plots", "figS3_test_sens.jpg"),
    p2,
    width = 5,
    height = 3.5,
    scale = 1,
    dpi = 300
)

readr::write_csv(plot_df,
                 here::here("output", "figS3_data.csv"))
