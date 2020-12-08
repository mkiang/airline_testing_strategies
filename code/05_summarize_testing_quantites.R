## 05_summarize_testing_quantites.R ----
##
## This collects important testing (and time-invariant) quantities from all
## the intermediate files and saves them into a single file.

## Imports ----
library(tidyverse)
library(here)
library(fs)
library(foreach)
library(doParallel)
source(here::here("code", "utils.R"))

testing_scenarios <- basename(fs::dir_ls(
    here::here("intermediate_files"),
    type = "directory",
    regexp = "pcr|rapid|no_testing|perfect"
))

## Calculate and summarize testing quantities ----
doParallel::registerDoParallel()
testing_results <- foreach::foreach(s = testing_scenarios) %dopar% {
    temp_x <-
        readRDS(here::here("data", sprintf("raw_simulations_%s.RDS", s)))
    
    ## Testing only occurs in a subset of time_steps so we don't need all of them
    x1 <- temp_x %>%
        dplyr::group_by(
            testing_type,
            testing_cat,
            prob_inf,
            prob_inf_cat,
            sens_type,
            symptom_screening,
            symptom_cat,
            risk_multiplier,
            risk_multi_cat,
            rapid_test_multiplier,
            rapid_test_cat,
            if_threshold,
            round,
            rep,
            sim_id
        ) %>%
        dplyr::summarize(
            n_test_false_pos_first = dplyr::case_when(
                ## No testing
                testing_type == "no_testing" ~ NA_integer_,
                ## Case of multiple tests
                testing_type == "rapid_same_day_pcr_three_days_after" ~
                    as.integer(n_test_false_pos[time_step == 70]),
                testing_type == "pcr_three_days_before_5_day_quarantine_pcr" ~
                    as.integer(n_test_false_pos[time_step == 70]),
                testing_type == "rapid_same_day_5_day_quarantine_pcr" ~
                    as.integer(n_test_false_pos[time_step == 70]),
                ## Case of just a single test
                TRUE ~ as.integer(max(n_test_false_pos, na.rm = TRUE))
            ),
            n_test_false_pos_second = dplyr::case_when(
                ## Case of multiple tests
                testing_type == "rapid_same_day_pcr_three_days_after" ~
                    as.integer(n_test_false_pos[time_step == 73]),
                testing_type == "pcr_three_days_before_5_day_quarantine_pcr" ~
                    as.integer(n_test_false_pos[time_step == 75]),
                testing_type == "rapid_same_day_5_day_quarantine_pcr" ~
                    as.integer(n_test_false_pos[time_step == 75]),
                ## Case of just a single test
                TRUE ~ NA_integer_
            ),
            n_test_true_pos_first = dplyr::case_when(
                ## No testing
                testing_type == "no_testing" ~ NA_integer_,
                ## Case of multiple tests
                testing_type == "rapid_same_day_pcr_three_days_after" ~
                    as.integer(n_test_true_pos[time_step == 70]),
                testing_type == "pcr_three_days_before_5_day_quarantine_pcr" ~
                    as.integer(n_test_true_pos[time_step == 70]),
                testing_type == "rapid_same_day_5_day_quarantine_pcr" ~
                    as.integer(n_test_true_pos[time_step == 70]),
                ## Case of just a single test
                TRUE ~ as.integer(max(n_test_true_pos, na.rm = TRUE))
            ),
            n_test_true_pos_second = dplyr::case_when(
                ## Case of multiple tests
                testing_type == "rapid_same_day_pcr_three_days_after" ~
                    as.integer(n_test_true_pos[time_step == 73]),
                testing_type == "pcr_three_days_before_5_day_quarantine_pcr" ~
                    as.integer(n_test_true_pos[time_step == 75]),
                testing_type == "rapid_same_day_5_day_quarantine_pcr" ~
                    as.integer(n_test_true_pos[time_step == 75]),
                ## Case of just a single test
                TRUE ~ NA_integer_
            ),
            n_test_false_pos = dplyr::case_when(
                ## No testing
                testing_type == "no_testing" ~ NA_integer_,
                ## Case of multiple tests
                testing_type == "rapid_same_day_pcr_three_days_after" ~
                    as.integer(n_test_false_pos[time_step == 70] +
                                   n_test_false_pos[time_step == 73]),
                testing_type == "pcr_three_days_before_5_day_quarantine_pcr" ~
                    as.integer(n_test_false_pos[time_step == 70] +
                                   n_test_false_pos[time_step == 75]),
                testing_type == "rapid_same_day_5_day_quarantine_pcr" ~
                    as.integer(n_test_false_pos[time_step == 70] +
                                   n_test_false_pos[time_step == 75]),
                ## Case of just a single test
                TRUE ~ as.integer(max(n_test_false_pos, na.rm = TRUE))
            ),
            n_test_true_pos = dplyr::case_when(
                ## No testing
                testing_type == "no_testing" ~ NA_integer_,
                ## Case of multiple tests
                testing_type == "rapid_same_day_pcr_three_days_after" ~
                    as.integer(n_test_true_pos[time_step == 70] +
                                   n_test_true_pos[time_step == 73]),
                testing_type == "pcr_three_days_before_5_day_quarantine_pcr" ~
                    as.integer(n_test_true_pos[time_step == 70] +
                                   n_test_true_pos[time_step == 75]),
                testing_type == "rapid_same_day_5_day_quarantine_pcr" ~
                    as.integer(n_test_true_pos[time_step == 70] +
                                   n_test_true_pos[time_step == 75]),
                ## Case of just a single test
                TRUE ~ as.integer(max(n_test_true_pos, na.rm = TRUE))
            ),
            n_infected_day_of_flight = n_infected_all[time_step == 70],
            n_active_infected_day_of_flight = n_active_infection[time_step == 70],
            n_active_infected_subclin_day_of_flight = n_active_infection_subclin[time_step == 70]
        ) %>%
        dplyr::ungroup() %>%
        dplyr::distinct() %>%
        dplyr::mutate(
            any_positive_test = n_test_false_pos + n_test_true_pos,
            ratio_false_true = n_test_false_pos / n_test_true_pos,
            ratio_true_false = n_test_true_pos / n_test_false_pos,
            ppv = n_test_true_pos / (n_test_true_pos + n_test_false_pos),
            frac_detected = n_test_true_pos / n_infected_day_of_flight,
            frac_active_detected = n_test_true_pos / n_active_infected_day_of_flight,
            time_step = 84
        )
    
    ## Calculate total infections
    x2 <- temp_x %>%
        dplyr::group_by(
            testing_type,
            testing_cat,
            prob_inf,
            prob_inf_cat,
            sens_type,
            symptom_screening,
            symptom_cat,
            risk_multiplier,
            risk_multi_cat,
            rapid_test_multiplier,
            rapid_test_cat,
            if_threshold,
            round,
            rep
        ) %>%
        dplyr::summarize(n_total_infections = n_susceptible[time_step == 84] -
                             n_susceptible[time_step == 67] +
                             n_infected_all[time_step == 67]) %>%
        dplyr::ungroup() %>%
        dplyr::distinct() %>%
        dplyr::mutate(time_step = 84)
    
    dplyr::full_join(x1, x2)
}
doParallel::stopImplicitCluster()
testing_results <- dplyr::bind_rows(testing_results)

saveRDS(testing_results,
        here::here("data", "all_testing_results.RDS"),
        compress = "xz")

## Summarize these results
testing_results <- testing_results %>%
    dplyr::group_by(
        testing_type,
        testing_cat,
        prob_inf,
        prob_inf_cat,
        sens_type,
        symptom_screening,
        symptom_cat,
        risk_multiplier,
        risk_multi_cat,
        rapid_test_multiplier,
        rapid_test_cat,
        if_threshold,
        time_step
    )

testing_summary <- dplyr::bind_rows(
    testing_results %>% summarize_results_column(any_positive_test),
    testing_results %>% summarize_results_column(n_test_true_pos),
    testing_results %>% summarize_results_column(n_test_false_pos),
    testing_results %>% summarize_results_column(n_test_true_pos_first),
    testing_results %>% summarize_results_column(n_test_false_pos_first),
    testing_results %>% summarize_results_column(n_test_true_pos_second),
    testing_results %>% summarize_results_column(n_test_false_pos_second),
    testing_results %>% summarize_results_column(n_active_infected_day_of_flight),
    testing_results %>% summarize_results_column(n_active_infected_subclin_day_of_flight),
    testing_results %>% summarize_results_column(n_infected_day_of_flight),
    testing_results %>% summarize_results_column(ratio_false_true),
    testing_results %>% summarize_results_column(ratio_true_false),
    testing_results %>% summarize_results_column(ppv),
    testing_results %>% summarize_results_column(frac_detected),
    testing_results %>% summarize_results_column(frac_active_detected),
    testing_results %>% summarize_results_column(n_total_infections)
) %>%
    dplyr::ungroup()

saveRDS(testing_summary,
        here::here("data", "summarized_testing_results.RDS"),
        compress = "xz")
