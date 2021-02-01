## 06_summarize_testing_quantites.R ----
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

## Different levels of adherence to symptom screening
adherence_level <- seq(0, 1, .2)

testing_dict <- list(
    no_testing = 70,
    pcr_two_days_before = 68,
    pcr_two_days_before_5_day_quarantine_pcr = c(70, 75),
    pcr_three_days_before = 68,
    pcr_three_days_before_5_day_quarantine_pcr = c(70, 75),
    pcr_three_days_before_7_day_quarantine_pcr = c(70, 77),
    pcr_three_days_before_14_day_quarantine_pcr = c(70, 84),
    pcr_five_days_before = 66, 
    pcr_five_days_before_5_day_quarantine_pcr = c(70, 75),
    pcr_seven_days_before = 64, 
    pcr_seven_days_before_5_day_quarantine_pcr = c(70, 75),
    rapid_test_same_day = 70,
    rapid_same_day_5_day_quarantine_pcr = c(70, 75),
    rapid_same_day_7_day_quarantine_pcr = c(70, 77),
    rapid_same_day_14_day_quarantine_pcr = c(70, 84),
    pcr_five_days_after = 75
)

## Calculate testing quantities ----
doParallel::registerDoParallel()
testing_results <- foreach::foreach(i = 1:NROW(testing_dict)) %dopar% {
    testing_type <- names(testing_dict)[i]
    time_steps <- testing_dict[[testing_type]]
    
    temp_x <- readRDS(here("data_raw", 
                           sprintf("raw_simulations_%s.RDS", 
                                   testing_type))) %>% 
        dplyr::group_by(
            testing_type,
            testing_cat,
            prob_inf,
            prob_inf_cat,
            sens_type,
            sens_cat,
            prop_subclin,
            prop_subclin_cat,
            symptom_screening,
            risk_multiplier,
            risk_multi_cat,
            rapid_test_multiplier,
            rapid_test_cat,
            if_threshold,
            round,
            rep,
            sim_id
        ) 
    
    ## Get day of flight infections and total infections
    temp_x_infs <- left_join(
        temp_x %>%
            filter(time_step == 70) %>%
            select(
                n_infected_day_of_flight = n_infected_all,
                n_active_infected_day_of_flight = n_active_infection,
                n_active_infected_subclin_day_of_flight = n_active_infection_subclin
            ),
        temp_x  %>%
            dplyr::summarize(n_total_infections = n_susceptible[time_step == 84] -
                                 n_susceptible[time_step == 67] +
                                 n_infected_all[time_step == 67]) %>%
            dplyr::ungroup() %>%
            dplyr::distinct()
    )
    
    ## Get testing statistics based on the scenario dictionary
    ## because some tests are *before* t-3 (when we started counting cumulative
    ## infections), we need to get the raw files instead of processed files. 
    temp_x_test <- collect_results(testing_type) %>% 
        dplyr::group_by(
            testing_type,
            prob_inf,
            sens_type,
            prop_subclin,
            symptom_screening,
            risk_multiplier,
            rapid_test_multiplier,
            if_threshold,
            round,
            rep
        ) %>% 
        filter(time_step %in% time_steps) 
    
    if (!has_name(temp_x_test, "n_test_false_pos")) {
        temp_x_test <- temp_x_test %>% 
            mutate(n_test_false_pos = NA_integer_,
                   n_test_true_pos = NA_integer_)
    }
    
    temp_x_test <- temp_x_test %>% 
        select(time_step, n_test_false_pos, n_test_true_pos) %>% 
        summarize(
            n_test_false_pos_first = case_when(
                testing_type == "no_testing"  ~ NA_integer_,
                n_distinct(time_steps) == 1 ~ as.integer(max(n_test_false_pos)),
                n_distinct(time_steps) == 2 ~ as.integer(n_test_false_pos[time_step == min(time_steps)])
            ),
            n_test_false_pos_second = case_when(
                testing_type == "no_testing"  ~ NA_integer_,
                n_distinct(time_steps) == 1 ~ NA_integer_,
                n_distinct(time_steps) == 2 ~ as.integer(n_test_false_pos[time_step == max(time_steps)])
            ),
            n_test_true_pos_first = case_when(
                testing_type == "no_testing"  ~ NA_integer_,
                n_distinct(time_steps) == 1 ~ as.integer(max(n_test_true_pos)),
                n_distinct(time_steps) == 2 ~ as.integer(n_test_true_pos[time_step == min(time_steps)])
            ),
            n_test_true_pos_second = case_when(
                testing_type == "no_testing"  ~ NA_integer_,
                n_distinct(time_steps) == 1 ~ NA_integer_,
                n_distinct(time_steps) == 2 ~ as.integer(n_test_true_pos[time_step == max(time_steps)])
            )
        ) %>%
        ungroup() %>% 
        distinct() %>%
        rowwise() %>% 
        mutate(
            n_test_false_pos = case_when(
                testing_type == "no_testing"  ~ NA_integer_,
                TRUE ~ sum(n_test_false_pos_first, n_test_false_pos_second, na.rm = TRUE)
            ),
            n_test_true_pos = case_when(
                testing_type == "no_testing"  ~ NA_integer_,
                TRUE ~ sum(n_test_true_pos_first, n_test_true_pos_second, na.rm = TRUE)
            )
        ) %>% 
        ungroup()
    
    rm(temp_x); gc2()
    
    ## Pivot wider to we can get "adherence" to symptom screening
    ## This shouldn't change anything but reviewers want it. 
    temp_x_wide <- temp_x_infs %>%
        left_join(temp_x_test) %>% 
        pivot_wider(
            id_cols = c(
                testing_type,
                testing_cat,
                prob_inf,
                prob_inf_cat,
                sens_type,
                sens_cat,
                prop_subclin,
                prop_subclin_cat,
                risk_multi_cat,
                risk_multiplier,
                rapid_test_multiplier,
                rapid_test_cat,
                if_threshold,
                round,
                rep,
                sim_id
            ),
            names_from = symptom_screening,
            values_from = n_infected_day_of_flight:n_test_true_pos
        )
    
    ## Loop through different levels of adherence to symptom screening. 
    temp_x_list <- vector("list", length = NROW(adherence_level))
    for (i in 1:NROW(adherence_level)) {
        a <- adherence_level[i]
        
        temp_x_list[[i]] <- temp_x_wide %>%
            transmute(
                testing_type,
                testing_cat,
                prob_inf,
                prob_inf_cat,
                sens_type,
                sens_cat,
                prop_subclin,
                prop_subclin_cat,
                risk_multiplier,
                risk_multi_cat,
                rapid_test_multiplier,
                rapid_test_cat,
                if_threshold,
                round,
                rep,
                sim_id,
                symptom_adherence = a,
                n_infected_day_of_flight = round(
                    (n_infected_day_of_flight_TRUE * a) +
                        (n_infected_day_of_flight_FALSE * (1 - a))
                ),
                n_active_infected_day_of_flight = round(
                    (n_active_infected_day_of_flight_TRUE * a) +
                        (n_active_infected_day_of_flight_FALSE * (1 - a))
                ),
                n_active_infected_subclin_day_of_flight = round(
                    (n_active_infected_subclin_day_of_flight_TRUE * a) +
                        (n_active_infected_subclin_day_of_flight_FALSE * (1 - a))
                ),
                n_total_infections = round((n_total_infections_TRUE * a) +
                                               (n_total_infections_FALSE * (1 - a))),
                
                n_test_false_pos_first = round(
                    (n_test_false_pos_first_TRUE * a) +
                        (n_test_false_pos_first_FALSE * (1 - a))
                ),
                n_test_false_pos_second = round(
                    (n_test_false_pos_second_TRUE * a) +
                        (n_test_false_pos_second_FALSE * (1 - a))
                ),
                n_test_true_pos_first = round(
                    (n_test_true_pos_first_TRUE * a) +
                        (n_test_true_pos_first_FALSE * (1 - a))
                ),
                n_test_true_pos_second = round(
                    (n_test_true_pos_second_TRUE * a) +
                        (n_test_true_pos_second_FALSE * (1 - a))
                ),
                n_test_false_pos = round((n_test_false_pos_TRUE * a) +
                                             (n_test_false_pos_FALSE * (1 - a))),
                n_test_true_pos = round((n_test_true_pos_TRUE * a) +
                                            (n_test_true_pos_FALSE * (1 - a)))
            )
    }
    
    ## Calculate testing quantities we are interested in but don't summarize
    temp_x_list %>%
        bind_rows() %>% 
        dplyr::mutate(
            any_positive_test = n_test_false_pos + n_test_true_pos,
            ratio_false_true = n_test_false_pos / n_test_true_pos,
            ratio_true_false = n_test_true_pos / n_test_false_pos,
            ppv = n_test_true_pos / (n_test_true_pos + n_test_false_pos),
            frac_detected = n_test_true_pos / n_infected_day_of_flight,
            frac_active_detected = n_test_true_pos / n_active_infected_day_of_flight,
            time_step = 84
        )
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
        sens_cat,
        prop_subclin,
        prop_subclin_cat,
        symptom_adherence,
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
