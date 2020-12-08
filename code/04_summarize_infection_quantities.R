## 04_summarize_infection_quantities.R ----
## 
## Once all the testing scenarios have run, we collect and summarize the 
## saved intermediate files into a single file. collect_and_munge_simulations()
## performs the heavy list here. Differences are calculated using the same
## simulation with identical random realizations (i.e., a testing scenario is
## compared to an identical non-testing scenario with the same symptomatic 
## draw, specificity draw, etc.).

## Imports ----
library(tidyverse)
library(here)
library(fs)
library(foreach)
library(doParallel)
source(here("code", "utils.R"))

testing_scenarios <- basename(
    fs::dir_ls(here("intermediate_files"),
               type = "directory",
               regexp = "pcr|rapid|no_testing|perfect")
    )

## Collect and save intermediate files ----
doParallel::registerDoParallel()
all_results <- foreach(s = testing_scenarios) %dopar% {
    temp_x <- collect_and_munge_simulations(s, results_x = "all")
    
    ### Save raw simulations ----
    saveRDS(temp_x, 
            here("data", sprintf("raw_simulations_%s.RDS", s)),
            compress = "xz")
}
doParallel::stopImplicitCluster()

## Get in the null results so we can calculate differences 
null_results <- readRDS(here("data", "raw_simulations_no_testing.RDS"))
no_symp <- null_results %>%
    filter(symptom_screening == FALSE) %>% 
    mutate(null_n_active_infection_clin = n_active_infection - n_active_infection_subclin) %>%
    select(
        prob_inf,
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

## Summarize across simulations ----
## Outcomes we want:
##      1. Cumulative infections over time
##      2. Cumulative infectious days over time
##      3. Number of new infections over time
##      4. Number of infections (weighted) over time
##      5. Number of infectious days (weighted) over time

## Calculate and summarize quantities of interest ----
doParallel::registerDoParallel()
summarized_results <- foreach(s = testing_scenarios) %dopar% {
    temp_x <- readRDS(here("data", sprintf("raw_simulations_%s.RDS", s)))
    
    ## Calculate number of active infections (clinical)
    temp_x <- temp_x %>% 
        mutate(n_active_infection_clin = n_active_infection - n_active_infection_subclin) 
    
    ## Join with null model (no testing, no symptom screening)
    temp_x <- temp_x %>% 
        left_join(no_symp)
    
    ## Absolute differences
    temp_x <- temp_x %>%
        mutate(
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
        mutate(
            rel_n_infected_all = (null_n_infected_all - n_infected_all) / null_n_infected_all,
            rel_n_infected_new = (null_n_infected_new - n_infected_new) / null_n_infected_new,
            rel_n_active_infection = (null_n_active_infection - n_active_infection) / null_n_active_infection,
            rel_n_active_infection_subclin = (null_n_active_infection_subclin - n_active_infection_subclin) / null_n_active_infection_subclin,
            rel_n_active_infection_clin = (null_n_active_infection_clin - n_active_infection_clin) / null_n_active_infection_clin,
            rel_n_infection_daily = (null_n_infection_daily - n_infection_daily) / null_n_infection_daily,
            rel_w_infection_daily = (null_w_infection_daily - w_infection_daily) / null_w_infection_daily,
            rel_cume_n_infected_new = (null_cume_n_infected_new - cume_n_infected_new) / null_cume_n_infected_new,
            rel_cume_n_infection_daily = (null_cume_n_infection_daily - cume_n_infection_daily) / null_cume_n_infection_daily,
            rel_cume_w_infection_daily = (null_cume_w_infection_daily - cume_w_infection_daily) / null_cume_w_infection_daily
        )
    
    ## Group these up and then use summarize function 
    temp_x <- temp_x %>% 
        group_by(
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
            time_step)
    
    ## summarize_results_column takes a single column and returns
    ## descriptive summary stats for that column (by the grouping variables)
    x <- bind_rows(
        temp_x %>% summarize_results_column(n_infected_all), 
        temp_x %>% summarize_results_column(n_infected_new), 
        temp_x %>% summarize_results_column(n_active_infection), 
        temp_x %>% summarize_results_column(n_active_infection_subclin), 
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
        ungroup()
}
doParallel::stopImplicitCluster()
summarized_results <- bind_rows(summarized_results)
saveRDS(summarized_results,
        here("data", "summarized_results.RDS"),
        compress = "xz")
