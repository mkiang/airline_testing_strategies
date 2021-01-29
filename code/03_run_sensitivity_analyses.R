## 03_run_sensitivity_analyses.R ----
## 
## Once the simulated populations have been created, we run all sensitivity
## analyses. See 02_run_main_simulations.R for more details as this script
## follows the same structure (differing only in the arguments passed to
## return_parameter_grid()). 

### Imports ----
library(tidyverse)
library(here)
library(fs)
library(config)
library(foreach)
library(doParallel)
source(here::here("code", "utils.R"))

### Constants / Unpack config file ----
cfig <- config::get(config = "dev")
n_cores <- cfig$n_cores
n_reps <- cfig$n_reps_per_round
subclin_infectious <- cfig$subclin_infectious
n_burnin <- cfig$n_burnin_days
day_of_flight <- cfig$day_of_flight
n_outcome_day <- cfig$n_outcome_day
days_quarantine <- cfig$days_quarantine

### Set up log folder
fs::dir_create(dir_logs())

### Create parameter grid ----
## These are sensitivity analyses that involve changing parameters of interest.
param_grid <- dplyr::bind_rows(
    return_parameter_grid(cfig = "sensitivity_risk_multiplier"),
    return_parameter_grid(cfig = "sensitivity_rt_multiplier"),
    return_parameter_grid(cfig = "sensitivity_sub_clin"), 
    return_parameter_grid(cfig = "sensitivity_test_sens"),
    return_sensitivity_parameter_grid()
    ) 

doParallel::registerDoParallel(cores = n_cores)
foreach::foreach(i = sample(1:NROW(param_grid))) %dopar% {
    ### Select parameters from the grid ----
    testing_type <- param_grid$testing_type[i]
    prob_inf <- param_grid$prob_inf[i]
    sens_type <- param_grid$sens_type[i]
    risk_multiplier <- param_grid$risk_multiplier[i]
    rapid_test_multiplier <- param_grid$rapid_test_multiplier[i]
    symptom_screening <- param_grid$symptom_screening[i]
    round <- param_grid$round[i]
    prop_subclin <- param_grid$prop_subclin[i]
    
    run_and_save_simulation(
        testing_type,
        prob_inf,
        sens_type,
        risk_multiplier,
        rapid_test_multiplier,
        symptom_screening,
        round,
        n_reps,
        prop_subclin,
        subclin_infectious,
        n_burnin,
        day_of_flight,
        n_outcome_day,
        days_quarantine
    )
}

## Close connections
doParallel::stopImplicitCluster()
