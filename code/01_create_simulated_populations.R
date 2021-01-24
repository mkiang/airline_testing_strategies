## 01_create_simulated_populations.R ----
## 
## This file creates the simulated cohorts and performs a burn-in period of
## ~2 months to ensure the simulation reaches a state of equilibrium before 
## the testing strategies are performed. All files are saved in the 
## intermediate_files folder along with the data necessary to reproduce the
## random state (e.g., .Random.seed objects). Comparisons between the testing
## strategy and a "no testing" scenario are *within* simulations using the 
## same random state and burned in population. 

## Imports ----
library(tidyverse)
library(here)
library(fs)
library(config)
library(foreach)
library(doParallel)
source(here::here("code", "utils.R"))

## Constants / Unpack config file ----
cfig <- config::get(config = "dev")
n_reps <- cfig$n_reps_per_round
n_rounds <- cfig$n_rounds_of_sims
n_pop <- cfig$n_passengers
prob_infs <- cfig$prob_inf / 1000000
prop_subclin <- cfig$prop_subclin
subclin_infectious <- cfig$subclin_infectious
n_cores <- cfig$n_cores
n_burnin <- cfig$n_burnin_days
days_quarantine <- cfig$days_quarantine
sens_type <- cfig$sens_type

## Create parameter grid of main analyses ----
param_grid <- expand.grid(
    round = 1:n_rounds,
    rep = 1:n_reps,
    prob_inf = prob_infs,
    prop_subclin = prop_subclin, 
    sens_type = sens_type, 
    stringsAsFactors = FALSE
)

## Bind the sensitivity analyses that involve new simulation states ----
param_grid <- bind_rows(
    param_grid,
    expand.grid(
        round = 1:config::get(config = "sensitivity_sub_clin")$n_rounds_of_sims,
        rep = 1:n_reps,
        prob_inf = prob_infs,
        prop_subclin = config::get(config = "sensitivity_sub_clin")$prop_subclin,
        sens_type = sens_type,
        stringsAsFactors = FALSE
    ),
    expand.grid(
        round = 1:config::get(config = "sensitivity_test_sens")$n_rounds_of_sims,
        rep = 1:n_reps,
        prob_inf = prob_infs,
        prop_subclin = prop_subclin,
        sens_type = config::get(config = "sensitivity_test_sens")$sens_type,
        stringsAsFactors = FALSE
    )
) %>% 
    filter(round %% 2 == 0)

## Remove files we've already created from the parameter grid
param_grid <- param_grid[with(param_grid,
                              !file.exists(
                                  return_sim_state_file(
                                      round,
                                      rep,
                                      prob_inf,
                                      prop_subclin,
                                      sens_type
                                  )
                              )), ]

doParallel::registerDoParallel(cores = n_cores)
foreach::foreach(i = sample(1:NROW(param_grid))) %dopar% {
    ### Select parameters from the grid ----
    round <- param_grid$round[i]
    rep <- param_grid$rep[i]
    prob_inf <- param_grid$prob_inf[i]
    sens_type <- param_grid$sens_type[i]
    prop_subclin <- param_grid$prop_subclin[i]
    
    ## Get the seed to save later
    initial_seed <- get_seed_alpha(list(Sys.time(), Sys.getpid()))
    
    ## Draw random variables (save later)
    set.seed(initial_seed)
    p_sens <- return_test_sensitivity(sens_type)
    p_spec <- draw_specificity(1)
    days_incubation <- draw_incubation(1)
    days_symptomatic <- draw_symptomatic(1)
    
    ## Draw infectiousness weights ----
    if_weights <- draw_infectiousness_weights(
        days_incubation = days_incubation,
        days_symptomatic = days_symptomatic
    )
    
    state_file_name <- return_sim_state_file(
        round, 
        rep, 
        prob_inf, 
        prop_subclin, 
        sens_type
        )
    
    fs::dir_create(dirname(state_file_name))
    
    if (!fs::file_exists(state_file_name)) {
        ## Create and burn in a simulated population ----
        ### Initialization ----
        sim_pop <- create_sim_pop(
            n_pop = n_pop,
            prop_subclin = prop_subclin,
            days_incubation = days_incubation,
            days_symptomatic = days_symptomatic
        )
        
        ### Burn-in ----
        sim_pop <- burnin_sim(
            sim_pop,
            n_iter = n_burnin - 1,
            if_weights = if_weights,
            days_quarantine = days_quarantine,
            prob_inf = prob_inf
        )
        
        ## Save everything we need to get back to current state
        x <- list(
            sim_pop = sim_pop,
            p_sens = p_sens,
            p_spec = p_spec,
            days_incubation = days_incubation,
            days_symptomatic = days_symptomatic,
            prop_subclin = prop_subclin, 
            sens_type = sens_type, 
            if_weights = if_weights,
            initial_seed = initial_seed,
            rseed = .Random.seed
        )
        
        saveRDS(x, state_file_name, compress = "xz")
    }
}
doParallel::stopImplicitCluster()
