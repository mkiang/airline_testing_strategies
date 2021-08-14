## 05_summarize_infection_quantities.R ----
##
## Once simulations have been collated and processed, we want our summary
## statistics aggregated over all simulations. Because of the way our simulation
## is set up, we can calculate "adherence" to symptom screening and quarantine
## as the weighted average between the non-screening/quarantine scenario and
## the quarantine scenario with perfect compliance. As a sensitivity, we
## take a bunch of different weighted averages to convey differing levels of
## adherence to both scenarios.

## Imports ----
library(tidyverse)
library(here)
library(fs)
library(foreach)
library(doParallel)
source(here::here("code", "utils.R"))

## Different levels of adherence to symptom screening
adherence_level <- seq(0, 1, .2)

## Get testing scenario names
testing_scenarios <- basename(fs::dir_ls(
    here::here("intermediate_files"),
    type = "directory",
    regexp = "pcr|rapid|no_testing|perfect"
))

## Get in the null results so we can calculate differences
null_results <-
    readRDS(here::here("data_raw", "raw_simulations_no_testing.RDS"))
no_symp <- null_results %>%
    dplyr::filter(symptom_screening == FALSE) %>%
    dplyr::mutate(null_n_active_infection_clin = n_active_infection - n_active_infection_subclin) %>%
    dplyr::select(
        prob_inf,
        prop_subclin,
        risk_multiplier,
        if_threshold,
        time_step,
        round,
        rep,
        null_n_infected_all = n_infected_all,
        null_n_infected_new = n_infected_new,
        null_n_infection_daily = n_infection_daily,
        null_w_infection_daily = w_infection_daily,
        null_n_active_infection = n_active_infection,
        null_n_active_infection_clin,
        null_n_active_infection_subclin = n_active_infection_subclin,
        null_cume_n_infection_daily = cume_n_infection_daily,
        null_cume_w_infection_daily = cume_w_infection_daily,
        null_cume_n_infected_new = cume_n_infected_new
    )

## Calculate outcomes along different levels of adherence ----
## Outcomes we want:
##      1. Cumulative infections over time
##      2. Cumulative infectious days over time
##      3. Number of new infections over time
##      4. Number of infections (weighted) over time
##      5. Number of infectious days (weighted) over time

doParallel::registerDoParallel()
temp_holder <- foreach::foreach(s = testing_scenarios) %dopar% {
    temp_x <-
        readRDS(here::here("data_raw", sprintf("raw_simulations_%s.RDS", s)))
    
    ## Calculate number of active infections (clinical)
    temp_x <- temp_x %>%
        dplyr::mutate(n_active_infection_clin = n_active_infection - n_active_infection_subclin)
    
    ## Join with null model (no testing, no symptom screening)
    temp_x <- temp_x %>%
        dplyr::left_join(no_symp)
    
    ## Absolute differences
    temp_x <- temp_x %>%
        dplyr::mutate(
            abs_n_infected_all = n_infected_all - null_n_infected_all,
            abs_n_infected_new = n_infected_new - null_n_infected_new,
            abs_n_active_infection = n_active_infection - null_n_active_infection,
            abs_n_active_infection_subclin = n_active_infection_subclin - null_n_active_infection_subclin,
            abs_n_active_infection_clin = n_active_infection_clin - null_n_active_infection_clin,
            abs_n_infection_daily = n_infection_daily - null_n_infection_daily,
            abs_w_infection_daily = w_infection_daily - null_w_infection_daily,
            abs_cume_n_infected_new = cume_n_infected_new - null_cume_n_infected_new,
            abs_cume_n_infection_daily = cume_n_infection_daily - null_cume_n_infection_daily,
            abs_cume_w_infection_daily = cume_w_infection_daily - null_cume_w_infection_daily
        )
    
    ## Relative differences
    temp_x <- temp_x %>%
        dplyr::mutate(
            rel_n_infected_all = (null_n_infected_all - n_infected_all) / null_n_infected_all,
            rel_n_infected_new = (null_n_infected_new - n_infected_new) / null_n_infected_new,
            rel_n_active_infection = (null_n_active_infection - n_active_infection) / null_n_active_infection,
            rel_n_active_infection_subclin = (
                null_n_active_infection_subclin - n_active_infection_subclin
            ) / null_n_active_infection_subclin,
            rel_n_active_infection_clin = (null_n_active_infection_clin - n_active_infection_clin) / null_n_active_infection_clin,
            rel_n_infection_daily = (null_n_infection_daily - n_infection_daily) / null_n_infection_daily,
            rel_w_infection_daily = (null_w_infection_daily - w_infection_daily) / null_w_infection_daily,
            rel_cume_n_infected_new = (null_cume_n_infected_new - cume_n_infected_new) / null_cume_n_infected_new,
            rel_cume_n_infection_daily = (null_cume_n_infection_daily - cume_n_infection_daily) / null_cume_n_infection_daily,
            rel_cume_w_infection_daily = (null_cume_w_infection_daily - cume_w_infection_daily) / null_cume_w_infection_daily
        )
    
    ## Pivot wider
    temp_x_wide <- temp_x %>%
        tidyr::pivot_wider(
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
                time_step,
                round,
                rep,
                sim_id
            ),
            names_from = symptom_screening,
            values_from = c(
                n_infected_all,
                n_infected_new,
                n_active_infection,
                n_active_infection_clin,
                n_infection_daily,
                w_infection_daily,
                cume_n_infected_new,
                cume_n_infection_daily,
                cume_w_infection_daily,
                abs_n_infected_all:rel_cume_w_infection_daily
            )
        )
    
    temp_x_list <- vector("list", length = NROW(adherence_level))
    for (i in 1:NROW(adherence_level)) {
        a <- adherence_level[i]
        
        temp_x_list[[i]] <- temp_x_wide %>%
            dplyr::transmute(
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
                time_step,
                round,
                rep,
                sim_id,
                symptom_adherence = a,
                n_infected_all = (n_infected_all_TRUE * a) + (n_infected_all_FALSE * (1 - a)),
                n_infected_new = (n_infected_new_TRUE * a) + (n_infected_new_FALSE * (1 - a)),
                n_active_infection = (n_active_infection_TRUE * a) + (n_active_infection_FALSE * (1 - a)),
                n_active_infection_clin = (n_active_infection_clin_TRUE * a) + (n_active_infection_clin_FALSE * (1 - a)),
                n_infection_daily = (n_infection_daily_TRUE * a) + (n_infection_daily_FALSE * (1 - a)),
                w_infection_daily = (w_infection_daily_TRUE * a) + (w_infection_daily_FALSE * (1 - a)),
                cume_n_infected_new = (cume_n_infected_new_TRUE * a) + (cume_n_infected_new_FALSE * (1 - a)),
                cume_n_infection_daily = (cume_n_infection_daily_TRUE * a) + (cume_n_infection_daily_FALSE * (1 - a)),
                cume_w_infection_daily = (cume_w_infection_daily_TRUE * a) + (cume_w_infection_daily_FALSE * (1 - a)),
                abs_n_infected_all = (abs_n_infected_all_TRUE * a) + (abs_n_infected_all_FALSE * (1 - a)),
                abs_n_infected_new = (abs_n_infected_new_TRUE * a) + (abs_n_infected_new_FALSE * (1 - a)),
                abs_n_active_infection = (abs_n_active_infection_TRUE * a) + (abs_n_active_infection_FALSE * (1 - a)),
                abs_n_active_infection_subclin = (abs_n_active_infection_subclin_TRUE * a) + (abs_n_active_infection_subclin_FALSE * (1 - a)),
                abs_n_active_infection_clin = (abs_n_active_infection_clin_TRUE * a) + (abs_n_active_infection_clin_FALSE * (1 - a)),
                abs_n_infection_daily = (abs_n_infection_daily_TRUE * a) + (abs_n_infection_daily_FALSE * (1 - a)),
                abs_w_infection_daily = (abs_w_infection_daily_TRUE * a) + (abs_w_infection_daily_FALSE * (1 - a)),
                abs_cume_n_infected_new = (abs_cume_n_infected_new_TRUE * a) + (abs_cume_n_infected_new_FALSE * (1 - a)),
                abs_cume_n_infection_daily = (abs_cume_n_infection_daily_TRUE * a) + (abs_cume_n_infection_daily_FALSE * (1 - a)),
                abs_cume_w_infection_daily = (abs_cume_w_infection_daily_TRUE * a) + (abs_cume_w_infection_daily_FALSE * (1 - a)),
                rel_n_infected_all = (rel_n_infected_all_TRUE * a) + (rel_n_infected_all_FALSE * (1 - a)),
                rel_n_infected_new = (rel_n_infected_new_TRUE * a) + (rel_n_infected_new_FALSE * (1 - a)),
                rel_n_active_infection = (rel_n_active_infection_TRUE * a) + (rel_n_active_infection_FALSE * (1 - a)),
                rel_n_active_infection_subclin = (rel_n_active_infection_subclin_TRUE * a) + (rel_n_active_infection_subclin_FALSE * (1 - a)),
                rel_n_active_infection_clin = (rel_n_active_infection_clin_TRUE * a) + (rel_n_active_infection_clin_FALSE * (1 - a)),
                rel_n_infection_daily = (rel_n_infection_daily_TRUE * a) + (rel_n_infection_daily_FALSE * (1 - a)),
                rel_w_infection_daily = (rel_w_infection_daily_TRUE * a) + (rel_w_infection_daily_FALSE * (1 - a)),
                rel_cume_n_infected_new = (rel_cume_n_infected_new_TRUE * a) + (rel_cume_n_infected_new_FALSE * (1 - a)),
                rel_cume_n_infection_daily = (rel_cume_n_infection_daily_TRUE * a) + (rel_cume_n_infection_daily_FALSE * (1 - a)),
                rel_cume_w_infection_daily = (rel_cume_w_infection_daily_TRUE * a) + (rel_cume_w_infection_daily_FALSE * (1 - a))
            )
    }
    
    saveRDS(dplyr::bind_rows(temp_x_list),
            here::here(
                "data_raw",
                sprintf("processed_simulations_wide_%s.RDS", s)
            ),
            compress = "xz")
}
doParallel::stopImplicitCluster()

## Summarize all results that do not incorporate quarantine ----
doParallel::registerDoParallel()
results_no_quarantine <-
    foreach::foreach(s = testing_scenarios[!grepl("quarantine", testing_scenarios)]) %dopar% {
        temp_x <-
            readRDS(here::here(
                "data_raw",
                sprintf("processed_simulations_wide_%s.RDS", s)
            ))
        
        ## Group these up and then use summarize function
        temp_x <- temp_x %>%
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
        
        ## summarize_results_column() takes a single column and returns
        ## descriptive summary stats for that column (by the grouping variables)
        dplyr::bind_rows(
            temp_x %>% summarize_results_column(n_infected_all),
            temp_x %>% summarize_results_column(n_infected_new),
            temp_x %>% summarize_results_column(n_active_infection),
            temp_x %>% summarize_results_column(n_active_infection_clin),
            temp_x %>% summarize_results_column(n_infection_daily),
            temp_x %>% summarize_results_column(w_infection_daily),
            temp_x %>% summarize_results_column(cume_n_infected_new),
            temp_x %>% summarize_results_column(cume_n_infection_daily),
            temp_x %>% summarize_results_column(cume_w_infection_daily),
            
            temp_x %>% summarize_results_column(abs_n_infected_all),
            temp_x %>% summarize_results_column(abs_n_infected_new),
            temp_x %>% summarize_results_column(abs_n_active_infection),
            temp_x %>% summarize_results_column(abs_n_active_infection_subclin),
            temp_x %>% summarize_results_column(abs_n_active_infection_clin),
            temp_x %>% summarize_results_column(abs_n_infection_daily),
            temp_x %>% summarize_results_column(abs_w_infection_daily),
            temp_x %>% summarize_results_column(abs_cume_n_infected_new),
            temp_x %>% summarize_results_column(abs_cume_n_infection_daily),
            temp_x %>% summarize_results_column(abs_cume_w_infection_daily),
            
            temp_x %>% summarize_results_column(rel_n_infected_all),
            temp_x %>% summarize_results_column(rel_n_infected_new),
            temp_x %>% summarize_results_column(rel_n_active_infection),
            temp_x %>% summarize_results_column(rel_n_active_infection_subclin),
            temp_x %>% summarize_results_column(rel_n_active_infection_clin),
            temp_x %>% summarize_results_column(rel_n_infection_daily),
            temp_x %>% summarize_results_column(rel_w_infection_daily),
            temp_x %>% summarize_results_column(rel_cume_n_infected_new),
            temp_x %>% summarize_results_column(rel_cume_n_infection_daily),
            temp_x %>% summarize_results_column(rel_cume_w_infection_daily)
        ) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(quarantine_adherence = 0)
    }
doParallel::stopImplicitCluster()

## Now like above, we want to take a weighted average using different weights
## to estimate different levels of adherence to quarantining.
quarantine_comparisons <- dplyr::bind_rows(
    expand.grid(
        base_case = "rapid_test_same_day",
        comparison_case = c(
            "rapid_same_day_5_day_quarantine_pcr",
            "rapid_same_day_7_day_quarantine_pcr",
            "rapid_same_day_14_day_quarantine_pcr"
        ),
        stringsAsFactors = FALSE
    ),
    expand.grid(
        base_case = "pcr_three_days_before",
        comparison_case = c(
            "pcr_three_days_before_5_day_quarantine_pcr",
            "pcr_three_days_before_7_day_quarantine_pcr",
            "pcr_three_days_before_14_day_quarantine_pcr"
        ),
        stringsAsFactors = FALSE
    ),
) %>%
    tibble::add_case(base_case = "pcr_five_days_before",
             comparison_case = "pcr_five_days_before_5_day_quarantine_pcr") %>%
    tibble::add_case(base_case = "pcr_seven_days_before",
             comparison_case = "pcr_seven_days_before_5_day_quarantine_pcr") %>%
    tibble::add_case(base_case = "pcr_two_days_before",
             comparison_case = "pcr_two_days_before_5_day_quarantine_pcr") %>%
    tibble::add_case(base_case = "pcr_five_days_after",
                     comparison_case = "5_day_quarantine_pcr_five_days_after")

doParallel::registerDoParallel()
results_quarantine <-
    foreach::foreach(i = 1:NROW(quarantine_comparisons)) %dopar% {
        b <- quarantine_comparisons$base_case[i]
        comp <- quarantine_comparisons$comparison_case[i]
        base_x <-
            readRDS(here::here(
                "data_raw",
                sprintf("processed_simulations_wide_%s.RDS", b)
            ))
        comp_x <-
            readRDS(here::here(
                "data_raw",
                sprintf("processed_simulations_wide_%s.RDS", comp)
            ))
        
        joined_x <- dplyr::left_join(
            comp_x,
            base_x %>%
                dplyr::select(-testing_type, -testing_cat),
            by = c(
                "prob_inf",
                "prob_inf_cat",
                "sens_type",
                "sens_cat",
                "prop_subclin",
                "prop_subclin_cat",
                "risk_multiplier",
                "risk_multi_cat",
                "rapid_test_multiplier",
                "rapid_test_cat",
                "if_threshold",
                "time_step",
                "round",
                "rep",
                "sim_id",
                "symptom_adherence"
            )
        )
        rm(base_x, comp_x)
        gc2()
        
        temp_x_list <- vector("list", length = NROW(adherence_level))
        for (i in 1:NROW(adherence_level)) {
            a <- adherence_level[i]
            
            temp_x_list[[i]] <- joined_x %>%
                dplyr::transmute(
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
                    time_step,
                    round,
                    rep,
                    sim_id,
                    symptom_adherence,
                    quarantine_adherence = a,
                    n_infected_all = (n_infected_all.x * a) + (n_infected_all.y * (1 - a)),
                    n_infected_new = (n_infected_new.x * a) + (n_infected_new.y * (1 - a)),
                    n_active_infection = (n_active_infection.x * a) + (n_active_infection.y * (1 - a)),
                    n_active_infection_clin = (n_active_infection_clin.x * a) + (n_active_infection_clin.y * (1 - a)),
                    n_infection_daily = (n_infection_daily.x * a) + (n_infection_daily.y * (1 - a)),
                    w_infection_daily = (w_infection_daily.x * a) + (w_infection_daily.y * (1 - a)),
                    cume_n_infected_new = (cume_n_infected_new.x * a) + (cume_n_infected_new.y * (1 - a)),
                    cume_n_infection_daily = (cume_n_infection_daily.x * a) + (cume_n_infection_daily.y * (1 - a)),
                    cume_w_infection_daily = (cume_w_infection_daily.x * a) + (cume_w_infection_daily.y * (1 - a)),
                    abs_n_infected_all = (abs_n_infected_all.x * a) + (abs_n_infected_all.y * (1 - a)),
                    abs_n_infected_new = (abs_n_infected_new.x * a) + (abs_n_infected_new.y * (1 - a)),
                    abs_n_active_infection = (abs_n_active_infection.x * a) + (abs_n_active_infection.y * (1 - a)),
                    abs_n_active_infection_subclin = (abs_n_active_infection_subclin.x * a) + (abs_n_active_infection_subclin.y * (1 - a)),
                    abs_n_active_infection_clin = (abs_n_active_infection_clin.x * a) + (abs_n_active_infection_clin.y * (1 - a)),
                    abs_n_infection_daily = (abs_n_infection_daily.x * a) + (abs_n_infection_daily.y * (1 - a)),
                    abs_w_infection_daily = (abs_w_infection_daily.x * a) + (abs_w_infection_daily.y * (1 - a)),
                    abs_cume_n_infected_new = (abs_cume_n_infected_new.x * a) + (abs_cume_n_infected_new.y * (1 - a)),
                    abs_cume_n_infection_daily = (abs_cume_n_infection_daily.x * a) + (abs_cume_n_infection_daily.y * (1 - a)),
                    abs_cume_w_infection_daily = (abs_cume_w_infection_daily.x * a) + (abs_cume_w_infection_daily.y * (1 - a)),
                    rel_n_infected_all = (rel_n_infected_all.x * a) + (rel_n_infected_all.y * (1 - a)),
                    rel_n_infected_new = (rel_n_infected_new.x * a) + (rel_n_infected_new.y * (1 - a)),
                    rel_n_active_infection = (rel_n_active_infection.x * a) + (rel_n_active_infection.y * (1 - a)),
                    rel_n_active_infection_subclin = (rel_n_active_infection_subclin.x * a) + (rel_n_active_infection_subclin.y * (1 - a)),
                    rel_n_active_infection_clin = (rel_n_active_infection_clin.x * a) + (rel_n_active_infection_clin.y * (1 - a)),
                    rel_n_infection_daily = (rel_n_infection_daily.x * a) + (rel_n_infection_daily.y * (1 - a)),
                    rel_w_infection_daily = (rel_w_infection_daily.x * a) + (rel_w_infection_daily.y * (1 - a)),
                    rel_cume_n_infected_new = (rel_cume_n_infected_new.x * a) + (rel_cume_n_infected_new.y * (1 - a)),
                    rel_cume_n_infection_daily = (rel_cume_n_infection_daily.x * a) + (rel_cume_n_infection_daily.y * (1 - a)),
                    rel_cume_w_infection_daily = (rel_cume_w_infection_daily.x * a) + (rel_cume_w_infection_daily.y * (1 - a))
                )
        }
        rm(joined_x)
        gc2()
        
        ## Group these up and then use summarize function
        temp_x <- temp_x_list %>%
            dplyr::bind_rows() %>%
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
                quarantine_adherence,
                risk_multiplier,
                risk_multi_cat,
                rapid_test_multiplier,
                rapid_test_cat,
                if_threshold,
                time_step
            )
        
        ## summarize_results_column() takes a single column and returns
        ## descriptive summary stats for that column (by the grouping variables)
        dplyr::bind_rows(
            temp_x %>% summarize_results_column(n_infected_all),
            temp_x %>% summarize_results_column(n_infected_new),
            temp_x %>% summarize_results_column(n_active_infection),
            temp_x %>% summarize_results_column(n_active_infection_clin),
            temp_x %>% summarize_results_column(n_infection_daily),
            temp_x %>% summarize_results_column(w_infection_daily),
            temp_x %>% summarize_results_column(cume_n_infected_new),
            temp_x %>% summarize_results_column(cume_n_infection_daily),
            temp_x %>% summarize_results_column(cume_w_infection_daily),
            
            temp_x %>% summarize_results_column(abs_n_infected_all),
            temp_x %>% summarize_results_column(abs_n_infected_new),
            temp_x %>% summarize_results_column(abs_n_active_infection),
            temp_x %>% summarize_results_column(abs_n_active_infection_subclin),
            temp_x %>% summarize_results_column(abs_n_active_infection_clin),
            temp_x %>% summarize_results_column(abs_n_infection_daily),
            temp_x %>% summarize_results_column(abs_w_infection_daily),
            temp_x %>% summarize_results_column(abs_cume_n_infected_new),
            temp_x %>% summarize_results_column(abs_cume_n_infection_daily),
            temp_x %>% summarize_results_column(abs_cume_w_infection_daily),
            
            temp_x %>% summarize_results_column(rel_n_infected_all),
            temp_x %>% summarize_results_column(rel_n_infected_new),
            temp_x %>% summarize_results_column(rel_n_active_infection),
            temp_x %>% summarize_results_column(rel_n_active_infection_subclin),
            temp_x %>% summarize_results_column(rel_n_active_infection_clin),
            temp_x %>% summarize_results_column(rel_n_infection_daily),
            temp_x %>% summarize_results_column(rel_w_infection_daily),
            temp_x %>% summarize_results_column(rel_cume_n_infected_new),
            temp_x %>% summarize_results_column(rel_cume_n_infection_daily),
            temp_x %>% summarize_results_column(rel_cume_w_infection_daily)
        ) %>%
            dplyr::ungroup()
    }
doParallel::stopImplicitCluster()

summarized_results <- dplyr::bind_rows(results_no_quarantine, results_quarantine)
saveRDS(summarized_results,
        here::here("data", "summarized_results.RDS"),
        compress = "xz")
