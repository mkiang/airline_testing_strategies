## Imports ----
library(tidyverse)
library(here)
library(fs)
source(here("code", "utils.R"))

## Constants ----
PROB_INF <- config::get()$primary_prob_inf / 1000000

## End of travel period infection stats
inf_df <- readRDS(here("data", "summarized_results.RDS")) %>%
    filter(
        time_step == 84, 
        if_threshold == 0,
        prob_inf == PROB_INF,
        risk_multiplier == 2,
        testing_type != "perfect_testing",
        rapid_test_multiplier == .9 |
            is.na(rapid_test_multiplier),
        sens_type != "median" | is.na(sens_type),
        prop_subclin == .3,
        symptom_adherence %in% c(0, .8),
        quarantine_adherence %in% c(0, .8)#,
        # testing_type %in% c(
        #     "no_testing",
        #     "pcr_three_days_before",
        #     "pcr_three_days_before_5_day_quarantine_pcr",
        #     "rapid_test_same_day",
        #     "rapid_same_day_5_day_quarantine_pcr",
        #     "pcr_five_days_after"
        # )
    ) %>%
    shift_time_steps() %>%
    categorize_metric() %>%
    filter(metric %in% c("cume_n_infection_daily",
                         "rel_cume_n_infection_daily"))

## Subset to the main analyses ----
inf_df <- bind_rows(
    ## No testing no symptom screening
    inf_df %>%
        filter(testing_type == "no_testing" &
                   symptom_adherence == 0),
    ## No quarantine 80% symptom screening adherence
    inf_df %>%
        filter(
            testing_type %in% c(
                "pcr_three_days_before",
                "rapid_test_same_day",
                "pcr_five_days_after",
                "pcr_five_days_before", 
                "pcr_seven_days_before",
                "pcr_two_days_before"
            ) & 
                symptom_adherence == .8
        ),
    ## With 80% quarantine and 80% symptom screening adherence
    inf_df %>% 
        filter(
            testing_type %in% c(
                "pcr_three_days_before_5_day_quarantine_pcr",
                "rapid_same_day_5_day_quarantine_pcr"
            ) & 
                symptom_adherence == .8 & 
                quarantine_adherence == .8
        )
) 

## Testing statistics 
test_df <- readRDS(here("data", "summarized_testing_results.RDS"))

test_df <- test_df %>%
    filter(
        risk_multiplier == 2,
        testing_type %in% c(
            "no_testing",
            "pcr_three_days_before",
            "pcr_three_days_before_5_day_quarantine_pcr",
            "rapid_test_same_day",
            "rapid_same_day_5_day_quarantine_pcr",
            "pcr_five_days_after",
            "pcr_two_days_before",
            "pcr_five_days_before",
            "pcr_seven_days_before"
        ),
        metric %in% c(
            "ratio_false_true",
            "any_positive_test",
            "n_test_true_pos",
            "n_test_false_pos",
            "n_active_infected_subclin_day_of_flight"
        ),
        rapid_test_multiplier == .9 | is.na(rapid_test_multiplier),
        prob_inf == PROB_INF,
        prop_subclin == .3,
        sens_type != "median" | is.na(sens_type),
        if_threshold == 0,
        symptom_adherence %in% c(0, .8)
    ) %>%
    categorize_metric()

test_df <- bind_rows(
    ## No testing no symptom screening
    test_df %>%
        filter(testing_type == "no_testing" &
                   symptom_adherence == 0),
    test_df %>%
        filter(testing_type != "no_testing" &
                   symptom_adherence == .8)
) 

## Day of flight statistics
dof_df <- readRDS(here("data", "summarized_results.RDS")) %>%
    filter(
        time_step == 70,
        if_threshold == 0,
        prob_inf == PROB_INF,
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
            "pcr_five_days_after",
            "pcr_two_days_before",
            "pcr_five_days_before",
            "pcr_seven_days_before"
        )
    ) %>%
    shift_time_steps() %>%
    categorize_metric() %>%
    filter(
        metric %in% c(
            "n_active_infection",
            "rel_n_active_infection",
            "n_active_infection_clin"
        )
    )

dof_df <- bind_rows(
    ## No testing no symptom screening
    dof_df %>%
        filter(testing_type == "no_testing" &
                   symptom_adherence == 0),
    ## No quarantine 80% symptom screening adherence
    dof_df %>%
        filter(
            testing_type %in% c(
                "pcr_three_days_before",
                "rapid_test_same_day",
                "pcr_five_days_after",
                "pcr_five_days_before", 
                "pcr_seven_days_before",
                "pcr_two_days_before"
            ) & 
                symptom_adherence == .8
        ),
    ## With 80% quarantine and 80% symptom screening adherence
    dof_df %>% 
        filter(
            testing_type %in% c(
                "pcr_three_days_before_5_day_quarantine_pcr",
                "rapid_same_day_5_day_quarantine_pcr"
            ) & 
                symptom_adherence == .8 & 
                quarantine_adherence == .8
        )
) 

## All main outcomes
all_df <- bind_rows(dof_df,
                    inf_df,
                    test_df) %>%
    arrange(testing_cat, metric_cat,) %>%
    mutate(print = case_when(
                  metric %in% c(
                      "cume_n_infection_daily",
                      "n_active_infection",
                      "n_active_infection_clin",
                      "n_active_infected_subclin_day_of_flight",
                      "any_positive_test",
                      "n_test_true_pos",
                      "n_test_false_pos"
                      
                  ) ~
                      sprintf("%i (95%% UI: %i, %i)",
                              round(mean),
                              round(p025),
                              round(p975)),
                  TRUE ~ sprintf(
                      "%0.2f (95%%: %0.2f, %0.2f)",
                      round(mean, 2),
                      round(p025, 2),
                      round(p975, 2)
                  )
              )) %>% 
    select(testing_cat, metric, metric_cat, print, everything())

write_csv(all_df, here("output", "99_main_manuscript_numbers.csv"))
saveRDS(all_df, here("output", "99_main_manuscript_numbers.RDS"))
