## 02_run_main_simulations.R ----
## 
## Once the simulated populations have been created, we run all scenarios 
## of interest and save the results in the intermediate files. The function
## return_parameter_grid() returns a list of parameters to sweep and checks
## for the existence of the intermediate file before running it (so redundant
## simulations are avoided). The function run_and_save_simulation() is a wrapper
## function that calls the correct testing scenario function. Details for
## both (and all other) functions are available in the utils.R file. 

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
param_grid <- return_parameter_grid(cfig = "dev")

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
closeAllConnections()
