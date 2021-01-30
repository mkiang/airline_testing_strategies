library(tidyverse)
library(here)
source(here("code", "utils.R"))
source(here("code", "mk_nytimes.R"))

plot_df <- sens_by_day_data() %>% 
    mutate(day = 1:n(),
           rt = .9 * upper)

p2 <- ggplot(plot_df, 
       aes(x = day, ymin = lower, ymax = upper, y = median)) + 
    geom_ribbon(alpha = .15) + 
    geom_point(alpha = .4) + 
    geom_line(data = plot_df, 
              aes(x = day, y = upper),
              color = "black",
              alpha = 1,
              size = .75) + 
    scale_x_continuous("Time from exposure, days", 
                       expand = c(.025, 0)) + 
    scale_y_continuous("Test sensitivity, %",
                       expand = c(.025, 0),
                       labels = function(x) round(x * 100)) + 
    mk_nytimes() 

ggsave(here("plots", "figS3_test_sens.pdf"), 
       p2,
       width = 5,
       height = 3.5,
       scale = 1,
       device = cairo_pdf)

ggsave(here("plots", "figS3_test_sens.jpg"), 
       p2,
       width = 5,
       height = 3.5,
       scale = 1,
       dpi = 300)

write_csv(plot_df,
          here("output", "figS3_data.csv"))
