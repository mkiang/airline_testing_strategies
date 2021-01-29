## Imports ----
library(truncnorm)
library(tidyverse)
library(data.table)
library(digest)
library(fs)
library(here)

## Comment out if you want "friendly" dplyr warnings.
options(dplyr.summarise.inform = FALSE)

## Testing simulations
## NOTE: Testing simulations are held in their own files because they are
## fairly long. See ./code/testing_scenarios for each file. There is a
## convenience wrapper called run_and_save_simulation() saved in the 
## ./code/testing_scenarios/run_and_save_simulation.R file. 
source(here("code", "testing_scenarios", "run_and_save_simulation.R"))

## Randomness helpers ----

#' Return a vector of n randomly generated non-symptomatic days
#'
#' For simplicity, this is really time spent asymptomatic, some of which is
#' incubation time (i.e., not infectious) and some of which is infectious 
#' (but again not symptomatic and therefore missed in symptomatic screening). 
#'
#' @param n size of randomly drawn incubation times
#' @param prop_pre_sx proportion of time pre-symptomatic but infectious?
#'
#' @return numeric vector
#' @example draw_symptomatic(1)
draw_incubation <- function(n, prop_pre_sx = 1/2) {
    ## This is a triangle distribution with median/mode at 3 (mean ~3.3) and
    ## bounds at [2, 5] after rounding.
    pre_sx_not_infectious <- round((rtri(n, 1.5, 5.49, 3)))
    
    ## This is a truncated normal with some proportion as pre-symptomatic
    pre_sx_infectious <- round(draw_duration_of_infectiousness(n) * prop_pre_sx)
    
    ## We return the sum of both days and keep the number of infectious pre-
    ## symptomatic days as the names
    res <- pre_sx_infectious + pre_sx_not_infectious
    names(res) <- pre_sx_infectious
    
    res
}

#' Return a vector for number of symptomatic days 
#' 
#' This returns a vector of symptomatic days, some of which will be infecitous
#' and some of which will be post-infectious but still symptomatic. 
#'
#' @param n 
#' @param days_test_pos_noninf 
#' @param prop_pre_sx 
draw_symptomatic <- function(n, days_test_pos_noninf = 14, prop_pre_sx = 1/2) {
    sx_infectious <- round(draw_duration_of_infectiousness(n) * prop_pre_sx)
    res <- sx_infectious + days_test_pos_noninf
    names(res) <- sx_infectious
    
    res
}

## Draw a vector for randomly generated duration of infectiousness
#' 
#' @param n 
#' @param mean 
#' @param a 
#' @param b 
#' @param sd 
draw_duration_of_infectiousness <- function(n, mean = 5, a = 1, b = Inf, sd = 1) {
    round(truncnorm::rtruncnorm(
        n = n,
        a = a,
        b = b,
        mean = mean,
        sd = sd
    ))
}

#' Return a trajectory of infectiousness for any given incubation/symptomatic
#'
#' Time-varying infectiousness weights. In our model, infections occur through
#' a daily static risk so there is no infectiousness curve. Below, we show
#' and example of how to generate an infectiousness curve in the comments. 
#'
#' @param days_incubation an int from array of draw_incubation()
#' @param days_symptomatic an int from array of draw_symptomatic()
#'
#' @return a named vector of infectiousness (of random length)
#' @example draw_infectiousness(draw_incubation(1), draw_symptomatic(1))
#' @export
# draw_infectiousness_weights <- function(days_incubation, days_symptomatic) {
#     ## Number of infectious days
#     mso <- truncnorm::rtruncnorm(1,
#                                  mean = 12.27248,
#                                  a = 9,
#                                  b = 15,
#                                  sd = 2)
#     d <- stats::dgamma(
#         seq(-(days_incubation - 1), days_symptomatic) + mso,
#         20.51651,
#         1.592124
#         )
#     d <- d / sum(d)
#     names(d) <- c(paste0("e", 1:days_incubation),
#                   paste0("l", 1:days_symptomatic))
#     
#     return(d)
# }
draw_infectiousness_weights <- function(days_incubation, days_symptomatic) {
    n_days_pre_sx_inf <- as.numeric(names(days_incubation))
    n_days_no_sx_no_inf <- days_incubation - n_days_pre_sx_inf
    n_days_sx_inf <- as.numeric(names(days_symptomatic))
    n_days_sx_no_inf <- days_symptomatic - n_days_sx_inf
    
    d <- c(
        rep(0, n_days_no_sx_no_inf),
        rep(1, n_days_pre_sx_inf),
        rep(1, n_days_sx_inf),
        rep(0, n_days_sx_no_inf)
    )
    
    names(d) <- c(paste0("e", 1:days_incubation),
                  paste0("l", 1:days_symptomatic))
    
    return(d)
}


#' Randomly draw specificity of the test, assumed to be between 99.5-100%%
#'
#' Mostly a convenience wrapper but keep it here in case we want to lower
#' specificity based on test type.
#'
#' @param n number of draws
#'
#' @return numeric between .98 and 1
#' @export
#'
#' @example draw_specificity()
draw_specificity <- function(n = 1,
                             lower = .995,
                             upper = 1,
                             mean = .998,
                             sd = .001) {
    truncnorm::rtruncnorm(
        n = n,
        a = lower,
        b = upper,
        mean = mean,
        sd = sd
    )
}

#' Return time-varying test sensitivity
#'
#' Returns the sensitivity of an RT-PCR by day of infection. Note, this can be
#' the median or upper/lower bounds of previous empirical estimates or it can
#' be perfect sensitivity each time ("perfect") or it can be random based on
#' empirical data.
#'
#' See sens_by_day_data() for notes about data source and correction.
#'
#' @param sens_type char in {"lower", "median", "upper", "perfect", "random"}
#' @return
#' @export
#'
#' @examples
return_test_sensitivity <- function(sens_type) {
    sens_by_day <- sens_by_day_data()
    
    if (is.na(sens_type)) {
        sens <- NA
    } else if (sens_type %in% names(sens_by_day)) {
        sens <- sens_by_day %>%
            dplyr::pull(sens_type)
    } else if (sens_type == "perfect") {
        sens <- rep(1, NROW(sens_by_day))
    } else if (sens_type == "random") {
        sens <- apply(sens_by_day, 1, function(x) {
            truncnorm::rtruncnorm(
                n = 1,
                a =  x[["lower"]],
                mean = x[["median"]],
                b = x[["upper"]],
                sd = (x[["upper"]] - x[["median"]]) / 2
            )
        })
    } else {
        stop(paste(
            c(
                "Not a valid sens_type. Must be one of:",
                names(sens_by_day),
                "perfect",
                "random"
            ),
            collapse = " "
        ))
    }
    
    sens
}


## Simulation helpers ----

#' Initialize a fully susceptible simulated population to follow
#'
#' Columns are:
#'     - ID: person ID
#'     - InfectionType: (If becomes infected) 0 = subclinical, 1 = clinical
#'     - State: 0 = Susceptible; 1 = Incubation; 2 = Symptomatic; 3 = Recovered
#'     - DaysIncubation: Days in incubation state
#'     - DaysSymptomatic: Days in symptomatic (infectious) state
#'     - NewInfection: Flag for infected on this day
#'     - DayofInfection: Day of initial infection
#'     - DayofDetection: Day of infection detection (PCR+ or symptom screen)
#'     - DayInObs: Day in OBSERVED state (not sure we need this)
#'
#' NOTE: May want to change this to all dplyr but for now keep data.table.
#'
#' @param n_pop number of people in this cohort
#' @param prop_subclin proportion of infections that are subclinical
#' @param days_incubation number of days of incubation (draw_incubation())
#' @param days_symptomatic number of days of symptomatic (draw_symptomatic())
#' @param as.dt return as a data.table (default = TRUE)
#'
#' @return a data.table (and tibble) at time 0 (i.e., all susceptible)
#' @export
#'
#' @example create_sim_pop(100, prop_subclin, draw_incubation(1), draw_symptomatic(1))
create_sim_pop <- function(n_pop,
                           prop_subclin,
                           days_incubation,
                           days_symptomatic) {
    res <- tibble::tibble(
        ID = as.character(1:n_pop),
        InfectionType = stats::rbinom(n_pop, 1, 1 - prop_subclin),
        State = 0,
        DaysIncubation = days_incubation,
        DaysPreSxInf = as.numeric(names(days_incubation)), 
        DaysSxInf = as.numeric(names(days_symptomatic)), 
        DaysSymptomatic = days_symptomatic,
        NewInfection = 0,
        DayOfInfection = 0,
        DayInObs = 0,
        DayOfDetection = 0
    )
    
    res <- data.table::as.data.table(res)
    data.table::setkey(res, ID)
    rownames(res) <- res$ID
    
    res
}

#' Takes relevant columns and time step and updates State column
#'
#' @param x sim_pop at time step t
#' @param t time step
#'
#' @return vector of State values {0, 1, 2, 3}
# update_state <- function(x, t) {
#     cols <- c("State",
#               "DayOfInfection",
#               "DaysIncubation",
#               "DaysSymptomatic")
#
#     x <- as.numeric(x[cols])
#     names(x) <- cols
#
#     if (x["State"] == 1 &
#         t >= (x["DaysIncubation"] + x["DayOfInfection"])) {
#         return(2) # Late Infectious period
#     }
#
#     if (x["State"] == 2 &
#         t >= (x["DaysIncubation"] +
#               x["DaysSymptomatic"] +
#               x["DayOfInfection"])) {
#         return(3) # Recovered
#     }
#
#     return(x["State"])
# }

#' Move the simulation forward one time step
#'
#' @param sim_pop the current iteration of the simulated population
#' @param t_step the current time step
#' @param if_weights infectious weights from draw_infectiousness_weights()
#' @param days_incubation single element from draw_incubation()
#' @param days_quarantine single element from draw_symptomatic()
#' @param p_community proportion of the community with active infection
#' @param prop_subclin proportion of infections that are subclinical
#' @param alpha_late_c proportion of symptomatic community that self-isolates
#' @param risk_multiplier risk multiplier for passengers
#' @param n_contacts number of daily contacts
#' @param beta_t infectiousness
#' @param subclin_infectious subclinical infectiousness discount
#' @param n_community number of people in the community
#'
#' @return
#' @export
increment_sim <- function(sim_pop,
                          t_step,
                          if_weights,
                          days_quarantine,
                          prob_inf) {
    ## Note: I changed the order of detected/quarantined below. I *think* the
    ## original order made it so DayInObs could never be 1 (i.e., it gets
    ## assigned 1 too early and then incremented at the end). Still need to
    ## test this order though.
    
    # Update quarantine day
    sim_pop[, DayInObs := ifelse(DayInObs <= days_quarantine &
                                     DayInObs > 0, DayInObs + 1, 0)]
    
    ## For newly detected cases, put them in quarantine
    sim_pop[DayOfDetection == t_step, DayInObs := 1]
    
    ## Reset NewInfection flag
    sim_pop[, NewInfection := 0]
    
    # S --> E
    # Get IDs of susceptible workers
    sWorkers <- sim_pop[State == 0 & DayInObs == 0, ID]
    
    ## Passengers have a static per-day risk of infection
    sim_pop[sWorkers, NewInfection := stats::rbinom(length(sWorkers), 1, prob_inf)]
    
    ## Update states
    sim_pop <- sim_pop %>%
        dplyr::mutate(State = dplyr::case_when(
            State == 1 & t_step >= (DaysIncubation + DayOfInfection) ~ 2,
            State == 2 &
                t_step >= (DaysIncubation + DaysSymptomatic + DayOfInfection) ~ 3,
            TRUE ~ State
        ))
    
    # Update infectious timelines
    sim_pop[NewInfection == 1, c("State", "DayOfInfection") := list(1, t_step + 1)]
    
    # Add time step for reference
    sim_pop$time_step <- t_step
    
    return(sim_pop)
}

#' Wrapper function to quickly increment a simulation (i.e., burn in)
#'
#' @param sim_pop simulated population
#' @param n_iter number of time steps to increment
#' @param ... parameters to pass to increment_sim()
#'
#' @return
#' @export
#'
#' @examples
burnin_sim <- function(sim_pop, n_iter = 50, ...) {
    for (i in 1:n_iter) {
        sim_pop <- increment_sim(sim_pop = sim_pop, t_step = i, ...)
    }
    return(sim_pop)
}

unpack_simulation_state <- function(f_name) {
    x <- readRDS(f_name)
    
    assign("sim_pop", data.table::copy(data.table::setDT(x$sim_pop)), envir = .GlobalEnv)
    assign("p_sens", x$p_sens, envir = .GlobalEnv)
    assign("p_spec", x$p_spec, envir = .GlobalEnv)
    assign("days_incubation", x$days_incubation, envir = .GlobalEnv)
    assign("days_symptomatic", x$days_symptomatic, envir = .GlobalEnv)
    assign("if_weights", x$if_weights, envir = .GlobalEnv)
    assign(".Random.seed", x$rseed, envir = .GlobalEnv)
}

## Parameter grid helpers ----
return_parameter_grid <- function(cfig = NA, new_only = TRUE) {
    cfig <- config::get(config = cfig)
    prob_infs <- cfig$prob_inf / 1000000
    n_rounds <- cfig$n_rounds_of_sims
    risk_multipliers <- cfig$risk_multiplier
    sens_type <- cfig$sens_type
    rapid_test_multipliers <- cfig$rt_sens_multiplier
    prop_subclin <- cfig$prop_subclin
    
    param_grid <- dplyr::bind_rows(
        ## Null model
        expand.grid(
            testing_type = "no_testing", 
            prob_inf = prob_infs,
            round = 1:n_rounds,
            sens_type = NA,
            risk_multiplier = risk_multipliers,
            rapid_test_multiplier = NA,
            symptom_screening = c(TRUE, FALSE),
            prop_subclin = prop_subclin, 
            stringsAsFactors = FALSE
        ),
        ## Null models
        expand.grid(
            testing_type = "perfect_testing", 
            prob_inf = prob_infs,
            round = 1:20,
            sens_type = NA,
            risk_multiplier = risk_multipliers,
            rapid_test_multiplier = NA,
            symptom_screening = c(TRUE, FALSE),
            prop_subclin = prop_subclin, 
            stringsAsFactors = FALSE
        ),
        ## Rapid tests
        expand.grid(
            testing_type = c(
                "rapid_test_same_day",
                "rapid_same_day_5_day_quarantine_pcr"
            ),
            prob_inf = prob_infs,
            sens_type = sens_type,
            risk_multiplier = risk_multipliers,
            rapid_test_multiplier = rapid_test_multipliers,
            symptom_screening = c(TRUE, FALSE),
            round = 1:n_rounds,
            prop_subclin = prop_subclin, 
            stringsAsFactors = FALSE
        ),
        ## PCR only
        expand.grid(
            testing_type = c(
                "pcr_five_days_after",
                "pcr_three_days_before",
                "pcr_three_days_before_5_day_quarantine_pcr"
            ),
            prob_inf = prob_infs,
            sens_type = sens_type,
            risk_multiplier = risk_multipliers,
            rapid_test_multiplier = NA,
            symptom_screening = c(TRUE, FALSE),
            round = 1:n_rounds,
            prop_subclin = prop_subclin, 
            stringsAsFactors = FALSE
        )
    )
    
    if (new_only) {
        param_grid <- param_grid[with(param_grid,
                                      !file.exists(
                                          return_sim_file_name(
                                              scenario_name = testing_type,
                                              symptom_screening = symptom_screening,
                                              prob_inf = prob_inf,
                                              prop_subclin = prop_subclin, 
                                              risk_multiplier = risk_multiplier,
                                              sens_type = sens_type, 
                                              rapid_test_multiplier = rapid_test_multiplier,
                                              round = round
                                          )
                                      )),]
    }
    param_grid
}

## Return sensitivity scenario parameter grid
return_sensitivity_parameter_grid <- function(new_only = TRUE) {
    cfig <- config::get(config = "sensitivity_scenarios")
    prob_infs <- cfig$prob_inf / 1000000
    n_rounds <- cfig$n_rounds_of_sims
    risk_multipliers <- cfig$risk_multiplier
    sens_type <- cfig$sens_type
    rapid_test_multipliers <- cfig$rt_sens_multiplier
    prop_subclin <- cfig$prop_subclin
    
    param_grid <- dplyr::bind_rows(
        ## PCR pre testing with different test timing (01a, 01b, 01c) or 
        ## quarantine duration (02a, 02b), or test timing + quarantine (02c/d/e)
        expand.grid(
            testing_type = c(
                "pcr_two_days_before",
                "pcr_five_days_before",
                "pcr_seven_days_before",
                "pcr_three_days_before_7_day_quarantine_pcr",
                "pcr_three_days_before_14_day_quarantine_pcr",
                "pcr_two_days_before_5_day_quarantine_pcr",
                "pcr_five_days_before_5_day_quarantine_pcr",
                "pcr_seven_days_before_5_day_quarantine_pcr"
            ),
            prob_inf = prob_infs,
            sens_type = sens_type,
            risk_multiplier = risk_multipliers,
            rapid_test_multiplier = NA,
            symptom_screening = c(TRUE, FALSE),
            round = 1:n_rounds,
            prop_subclin = prop_subclin, 
            stringsAsFactors = FALSE
        ),
        
        ## Rapid tests with different quarantine periods (04a, 04b)
        expand.grid(
            testing_type = c(
                "rapid_same_day_7_day_quarantine_pcr",
                "rapid_same_day_14_day_quarantine_pcr"
            ),
            prob_inf = prob_infs,
            sens_type = sens_type,
            risk_multiplier = risk_multipliers,
            rapid_test_multiplier = rapid_test_multipliers,
            symptom_screening = c(TRUE, FALSE),
            round = 1:n_rounds,
            prop_subclin = prop_subclin, 
            stringsAsFactors = FALSE
        )
    )
    
    if (new_only) {
        param_grid <- param_grid[with(param_grid,
                                      !file.exists(
                                          return_sim_file_name(
                                              scenario_name = testing_type,
                                              symptom_screening = symptom_screening,
                                              prob_inf = prob_inf,
                                              prop_subclin = prop_subclin, 
                                              risk_multiplier = risk_multiplier,
                                              sens_type = sens_type, 
                                              rapid_test_multiplier = rapid_test_multiplier,
                                              round = round
                                          )
                                      )),]
    }
    param_grid
}


## Testing helpers ----

#' Tested infected population (true positives)
#'
#' @param sim_pop_pos subset of sim_pop that are incubaton or symptomatic
#' @param t time step
#' @param p_sens vector of day-specific test sensitivity
#' @param days_incubation days of incubation
#' @param days_symptomatic days of symptomatic
#' @param n_delay test reporting delays
#' @param min_days days before symptoms where test can be positive
#'
#' @return subset of sim_pop
test_pos <- function(sim_pop_pos,
                     t,
                     p_sens,
                     days_incubation,
                     days_symptomatic,
                     n_delay = 0,
                     min_days = 5) {
    p_i <-
        ceiling(apply(sim_pop_pos, 1, function(y)
            day_in_infection(y, t)))
    
    # can only test positive X days prior to symptom onset
    toffset <- days_incubation - min_days
    p_i <- p_i - toffset
    
    sTested_tp <-
        sapply(p_i, function(y)
            ifelse(y > 0 & y <= length(p_sens),
                   stats::rbinom(1, 1, p_sens[y]),
                   0))
    
    sim_pop_pos[sTested_tp == 1, DayOfDetection := t + n_delay]
    pos_ids <- sim_pop_pos[sTested_tp == 1, ID]
    
    return(list(sim_pop_pos, pos_ids))
}

#' Test noninfected population (negative positives)
#'
#' @param sim_pop_neg subset of sim_pop of only uninfected
#' @param t time step
#' @param p_spec numeric of test specificity
#' @param n_delay test reporting delay
#'
#' @return subset of sim_pop
test_neg <- function(sim_pop_neg,
                     t,
                     p_spec,
                     n_delay = 0) {
    sTested_fp <- stats::rbinom(nrow(sim_pop_neg), 1, 1 - p_spec)
    sim_pop_neg[sTested_fp == 1, DayOfDetection := t + n_delay]
    neg_ids <- sim_pop_neg[sTested_fp == 1, ID]
    
    return(list(sim_pop_neg, neg_ids))
}



#' Test the entire simulated population
#'
#' HACK: Note that this incorporates truly horrific hacks because data.table
#' in-place modification doesn't like it when you add new columns. As a result
#' I'm returning lists for everything. May swap to dplyr at some point.
#'
#' @param sim_pop simulated population
#' @param t time step
#' @param p_sens vector of day-specific test sensitivity
#' @param p_spec numeric of test specificity
#' @param days_incubation number of days of incubation
#' @param days_symptomatic number of days of symptomatic
#' @param n_delay reporting delays
#' @param min_days minimum number of days before detectable virus
#'
#' @return
#' @export
test_pop <- function(sim_pop,
                     t,
                     p_sens,
                     p_spec,
                     days_incubation,
                     days_symptomatic,
                     n_delay = 0,
                     min_days = 5) {
    ## Test the infected
    ## Subset to infected IDs
    pids <-
        sim_pop[State %in% 1:2 & DayOfDetection < DayOfInfection, ID]
    
    test_pos <- test_pos(
        sim_pop[pids,],
        t,
        p_sens = p_sens,
        days_incubation = days_incubation,
        days_symptomatic = days_symptomatic,
        n_delay = n_delay,
        min_days = min_days
    )
    
    ## Replace infected with their rest results and save true positives
    sim_pop[pids,] <- test_pos[[1]]
    
    ## Test susceptible population
    test_neg <- test_neg(sim_pop[State == 0,],
                         t,
                         p_spec = p_spec,
                         n_delay = 0)
    sim_pop[State == 0,] <- test_neg[[1]]
    
    ## Test recovered
    rids <-
        sim_pop[State == 3 & DayOfDetection < DayOfInfection, ID]
    test_neg_recovered <- test_neg(sim_pop[rids,],
                                   t,
                                   p_spec = p_spec,
                                   n_delay = 0)
    sim_pop[rids,] <- test_neg_recovered[[1]]
    
    return(list(
        sim_pop = sim_pop,
        sTested_tp = test_pos[[2]],
        sTested_fp = unique(c(test_neg[[2]], test_neg_recovered[[2]]))
    ))
}

rapid_antigen_testing <- function(sim_pop,
                                  t,
                                  p_sens,
                                  p_spec,
                                  rapid_antigen_multiplier,
                                  days_incubation,
                                  days_symptomatic,
                                  n_delay = 0,
                                  min_days = 5) {
    x <- data.table::copy(sim_pop)
    
    ## Remove everybody who is in quarantine
    wids <- x[DayInObs == 0 & DayOfDetection < t, ID]
    
    ## Test the rest
    tested_pop <- test_pop(
        x[wids,],
        t = t,
        p_sens = p_sens * rapid_antigen_multiplier,
        p_spec = p_spec,
        days_incubation = days_incubation,
        days_symptomatic = days_symptomatic,
        n_delay = n_delay,
        min_days = min_days
    )
    
    ## Put new results back
    x[wids,] <- tested_pop[["sim_pop"]]
    
    ## Add columns for the false/true positive flags
    x %>%
        dplyr::mutate(
            sTested_tp = ifelse(ID %in% tested_pop[["sTested_tp"]], 1, 0),
            sTested_fp = ifelse(ID %in% tested_pop[["sTested_fp"]], 1, 0)
        )
}

pcr_testing <- function(sim_pop,
                        t,
                        p_sens,
                        p_spec,
                        days_incubation,
                        days_symptomatic,
                        n_delay = 0,
                        min_days = 5) {
    ## The test_pop function doesn't work if we introduce new columns so
    ## we need to check for previous testing columns.
    x <- data.table::copy(sim_pop)
    
    if (tibble::has_name(x, "sTested_tp") |
        tibble::has_name(x, "sTested_fp")) {
        if (any(!is.na(c(x$sTested_tp, x$sTested_fp)))) {
            stop("This population has already been tested before.")
        } else {
            stop("Remove sTested_tp and sTested_fp.")
        }
    }
    
    ## Remove everybody who is in quarantine
    wids <- x[DayInObs == 0 & DayOfDetection < t, ID]
    
    ## Test the rest
    tested_pop <- test_pop(
        x[wids,],
        t = t,
        p_sens = p_sens,
        p_spec = p_spec,
        days_incubation = days_incubation,
        days_symptomatic = days_symptomatic,
        n_delay = n_delay,
        min_days = min_days
    )
    
    ## Put new results back
    x[wids,] <- tested_pop[["sim_pop"]]
    
    ## Add columns for the false/true positive flags
    x %>%
        dplyr::mutate(
            sTested_tp = ifelse(ID %in% tested_pop[["sTested_tp"]], 1, 0),
            sTested_fp = ifelse(ID %in% tested_pop[["sTested_fp"]], 1, 0)
        )
}

#' Wrapper for disassembling the sim_pop, testing a subpop, then reassembling
#' the dataframe again.
#'
#' TODO: Right now, data.table is giving me issues with
#' in-place mods so I convert sim_pop to a normal df,
#' split up the population into to-be-tested and not-tested,
#' then convert sim_pop_a back to a dt so I can test the
#' sim_pop_a population, then put the dataset back together
#' and convert it back to a dt.
#'
#' This seems unnecessarily convoluted but fixing it would
#' require refactoring all the testing functions from
#' data.table to dplyr (I think).
#'
#' @param sim_pop
#' @param k
#' @param n_subsets
#' @param subset_to_test
#' @param p_sens
#' @param p_spec
#' @param days_incubation
#' @param days_symptomatic
#' @param n_delay
#' @param min_days
#'
#' @return
#' @export
#'
#' @examples
pcr_test_subset_simpop <- function(sim_pop,
                                   k,
                                   n_subsets = 2,
                                   subset_to_test = 1,
                                   p_sens,
                                   p_spec,
                                   days_incubation,
                                   days_symptomatic,
                                   n_delay = 0,
                                   min_days = 5) {
    test_ids <-
        return_portion_of_ids(sim_pop, n_subsets, subset_to_test)
    
    data.table::setDF(sim_pop)
    sim_pop_to_test <- sim_pop %>%
        dplyr::filter(ID %in% test_ids)
    sim_pop_no_test <- sim_pop %>%
        dplyr::filter(!(ID %in% test_ids))
    
    ## NOTE: YOU HAVE TO REMOVE THE EXTRA COLUMNS OR DATA.TABLE
    ## WON'T LET YOU RUN THE TESTING CODE.
    if (tibble::has_name(sim_pop_to_test, "sTested_tp") |
        tibble::has_name(sim_pop_to_test, "sTested_fp")) {
        if (any(!is.na(
            c(
                sim_pop_to_test$sTested_tp,
                sim_pop_to_test$sTested_fp
            )
        ))) {
            stop("This population has already been tested before.")
        } else {
            sim_pop_to_test <- sim_pop_to_test %>%
                dplyr::select(-sTested_tp,
                              -sTested_fp)
        }
    }
    
    ## Test the ones we need to test
    sim_pop_to_test <- set_id_as_key(sim_pop_to_test)
    sim_pop_to_test <- pcr_testing(sim_pop_to_test,
                                   k,
                                   p_sens,
                                   p_spec,
                                   days_incubation,
                                   days_symptomatic)
    data.table::setDF(sim_pop_to_test)
    
    ## Put it back together again and reset the key
    sim_pop <- dplyr::bind_rows(sim_pop_to_test,
                                sim_pop_no_test)
    sim_pop <- set_id_as_key(sim_pop)
    
    sim_pop
}

clear_old_tests <- function(sim_pop) {
    sim_pop[, c("sTested_fp", "sTested_tp") := NULL]
}

## Infections helpers ----
calculate_daily_infectiousness <- function(sim_pop,
                                           if_weights,
                                           subclin_infectious,
                                           fast_mode = FALSE,
                                           min_if_weight = 0,
                                           keep_cols = FALSE) {
    ## Not an infectious day unless meets minimum threshold
    if_binary <- (if_weights > min_if_weight) + 0
    
    ## Reweigh infectiousness limited to just infectious days (at threshold)
    if_weights_x <- ifelse(if_weights <= min_if_weight, 0, if_weights)
    if_weights_x <- if_weights_x / sum(if_weights_x)
    
    ## Prep the dataframe by adding number where they are in infectious
    ## trajectory and number of symptomatic days
    x <- sim_pop %>%
        add_day_since_infection()
    
    if (fast_mode) {
        x <- x %>%
            dplyr::filter(infection_day > 0, State != 3)
    }
    
    ## On this day, are they infectious (0, .5, 1) and if so, how infectious?
    x <- x %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            n_infection_daily = ifelse(
                infection_day > 0 &
                    infection_day <= (DaysIncubation + DaysSymptomatic) &
                    DayInObs == 0 &
                    !is.na(infection_day),
                if_binary[infection_day],
                0
            ),
            w_infection_daily = ifelse(
                infection_day > 0 &
                    infection_day <= (DaysIncubation + DaysSymptomatic) &
                    DayInObs == 0 &
                    !is.na(infection_day),
                if_weights_x[infection_day],
                0
            )
        ) %>%
        dplyr::ungroup() %>%
        ## Reduce subclinical infectiousness
        dplyr::mutate(
            w_infection_daily = ifelse(
                InfectionType == 0,
                subclin_infectious * w_infection_daily,
                w_infection_daily
            )
        ) %>%
        ## Reduce the number of infectious days for subclinical
        dplyr::mutate(
            n_infection_daily = ifelse(
                InfectionType == 0,
                subclin_infectious * n_infection_daily,
                n_infection_daily
            )
        )
    
    if (keep_cols) {
        x
    }
    
    x %>% dplyr::select(-infection_day)
}

## Data manipulation helpers ----

#' Add day since initial infection column
#'
#' @param sim_pop data frame with simulated population
#'
#' @return same data frame with a new column (infection_day)
add_day_since_infection <-  function(sim_pop) {
    ## On this time step, where in the infectious trajectory are they?
    sim_pop %>%
        dplyr::mutate(infection_day = time_step - ifelse(DayOfInfection > 0,
                                                         DayOfInfection,
                                                         NA))
}

#' Add number of symptomatic days to dataframe as a new column
#'
#' @param sim_pop
#'
#' @return same data frame with a new column (n_symptomatic_days)
# add_symptomatic_days <- function(sim_pop, day_of_flight) {
#     sim_pop %>%
#         mutate(
#             n_symptomatic_days = case_when(
#                 ## If you were quarantined, you only contribute as many days as
#                 ## your post-flight, non-quarantined days
#                 DayOfDetection > 0 ~ max(DayOfDetection - day_of_flight, 0),
#
#                 ## Remove infection days from those that would have been
#                 ## recovered by the day of the flight
#                 # (day_of_flight - DaysIncubation - DaysSymptomatic) >=
#                 #     DayOfInfection ~ 0,
#
#                 ## Current infection day cannot be after you would have recovered
#                 infection_day > 0 &
#                     infection_day <= (DaysIncubation + DaysSymptomatic) ~
#                     infection_day,
#
#                 ## Else you weren't infected
#                 TRUE ~ 0
#             )
#         )
# }
add_symptomatic_days <- function(sim_pop) {
    sim_pop %>%
        dplyr::mutate(
            n_symptomatic_days = dplyr::case_when(
                ## If they're positive, then number of days up until test
                DayOfInfection > 0 & DayOfDetection > 0 ~
                    DayOfDetection - DayOfInfection,
                
                ## Don't count those that are recovered
                infection_day > 0 &
                    infection_day <= (DaysIncubation + DaysSymptomatic) ~
                    infection_day,
                
                ## Else you weren't infected
                TRUE ~ 0
            )
        )
}


#' Add State_combined column
#'
#' A convenience wrapper that keeps track of subclinical and clinical states
#'
#'    0.0 - Susceptible
#'    1.0 - Subclinical incubation
#'    1.1 - Clinical incubation
#'    2.0 - Subclinical symptomatic
#'    2.1 - Clinical symptomatic
#'    3.0 - Recovered
#'
#' @param sim_pop Simulated population from create_sim_pop()
#'
#' @return sim_pop data frame with additional State_combined column
#' @export
add_State_combined <- function(sim_pop) {
    sim_pop %>%
        dplyr::mutate(State_combined = ifelse(State %in% 1:2,
                                              as.numeric(
                                                  sprintf("%i.%i", State, InfectionType)
                                              ),
                                              State))
}

add_postflight_inf_days <- function(sim_pop, day_of_flight) {
    sim_pop %>%
        dplyr::mutate(postflight_inf_days =
                          DaysSymptomatic + DaysIncubation + DayOfInfection -
                          day_of_flight) %>%
        dplyr::mutate(postflight_inf_days = ifelse(postflight_inf_days < 0,
                                                   0,
                                                   postflight_inf_days))
}

subset_to_inf_passengers <- function(sim_pop) {
    sim_pop %>%
        dplyr::filter(
            DayOfInfection >= (n_burnin - days_incubation -
                                   days_symptomatic + 1),
            DayOfInfection <= 120
        )
}

## Misc. helpers ----
return_sim_file_name <- function(scenario_name,
                                 symptom_screening,
                                 prob_inf,
                                 prop_subclin, 
                                 risk_multiplier,
                                 rapid_test_multiplier,
                                 sens_type, 
                                 round) {
    sprintf(
        "%s/%s/%03d-%s-prob_inf_%03d-asx_%02d-risk_multi_%i-sens_type_%s-rt_multi_%i-results.RDS",
        dir_results(scenario_name),
        ifelse(symptom_screening, "with_screening", "no_screening"),
        round,
        scenario_name,
        prob_inf * 1000000,
        prop_subclin * 100, 
        risk_multiplier,
        sens_type, 
        rapid_test_multiplier * 100
    )
}

return_sim_state_file <- function(round, rep, prob_inf, prop_subclin) {
    sprintf(
        "%s/round_%02d/round_%02d-rep_%03d-prob_inf_%03d-asx_%02d-post_burnin_simulation_state.RDS",
        dir_results("simulation_states"),
        round,
        round,
        rep,
        prob_inf * 1000000,
        prop_subclin * 100
    )
}

print_sim_start <- function(scenario_name,
                            symptom_screening,
                            prob_inf,
                            sens_type,
                            risk_multiplier,
                            rapid_test_multiplier,
                            prop_subclin, 
                            round) {
    print(
        sprintf(
            paste0(
                "Starting (%s): ",
                "test: %s, ",
                "round: %s, ",
                "sens: %s, ",
                "symptom_testing: %s, ",
                "prob_inf: %s, ",
                "prop_subclin: %s, ", 
                "risk_multi: %s, ",
                "rt_multi: %s"
            ),
            Sys.time(),
            scenario_name,
            round,
            sens_type,
            symptom_screening,
            prob_inf,
            prop_subclin, 
            risk_multiplier,
            rapid_test_multiplier
        )
    )
}

print_sim_rep <- function(rep,
                          days_incubation,
                          days_symptomatic,
                          prob_inf,
                          prop_subclin) {
    print(
        sprintf(
            "Rep %s (%s): days_inc: %s, days_sympt: %s, prob_inf: %s, prop_subclin: %s",
            rep,
            Sys.time(),
            days_incubation,
            days_symptomatic,
            prob_inf, 
            prop_subclin
        )
    )
}


set_id_as_key <- function(res) {
    res <- data.table::as.data.table(res)
    data.table::setkey(res, ID)
    rownames(res) <- res$ID
    res
}

#' For splitting up simulated pop into chunks systematically
#'
#' Given a sim_pop, will split up IDs into specified number of approximately-
#' equal rows and return one of the chunks. This is useful when we want
#' to systematically split up the population so we can test different subsets
#' on different days without overlap.
#'
#' @param sim_pop a simulated population
#' @param split_into_n_pieces int of number chunks to return
#' @param return_piece int of the chunk to return
#'
#' @return a vector of IDs
#' @export
#'
#' @examples
#' \dontrun{
#' identical(
#' sort(sim_pop$ID),
#' sort(c(
#'     return_portion_of_ids(sim_pop, n_split = 3, return_chunk = 1),
#'     return_portion_of_ids(sim_pop, n_split = 3, return_chunk = 2),
#'     return_portion_of_ids(sim_pop, n_split = 3, return_chunk = 3)
#' ))
#' )
#' # TRUE
#' }
return_portion_of_ids <- function(sim_pop,
                                  n_split = 2,
                                  return_chunk = 1) {
    x <- sim_pop$ID
    
    split(x, sort(1:length(x) %% n_split))[[return_chunk]]
}

dir_results <- function(scenario_name) {
    here::here("intermediate_files",
         # "simulations",
         scenario_name) #,
    # "sim_results")
}

# dir_sim_object <- function(scenario_name) {
#     here("intermediate_files",
#          "simulations",
#          scenario_name,
#          "sim_objects")
# }

dir_logs <- function() {
    here::here("intermediate_files", "logs")
}

#' Wrapper to return current day of infection
#'
#' @param x a single row from sim_pop (i.e., meant to be used in an apply())
#' @param t time step
#'
#' @return numeric
day_in_infection <- function(x, t) {
    cols <- c("State", "DayOfInfection")
    x <- as.numeric(x[cols])
    names(x) <- cols
    
    if (x["State"] %in% c(1, 2)) {
        return(t - x["DayOfInfection"] + 1)
    }
    
    return(NA)
}

#' Given any object, return an int
#'
#' Taken from https://github.com/etchin/covid-testing/blob/master/microsims.R
#'
#' @param x any valid R object (including NA, NULL, etc)
#'
#' @return int
#' @export
get_seed_alpha <- function(x) {
    hexval <- paste0("0x", digest::digest(x, "crc32"))
    intval <- utils::type.convert(hexval) %% .Machine$integer.max
    return(intval)
}


#' Return the corrected test sensitivity by day data.
#'
#' NOTE: This should be a dot function if we turn this into a package.
#'
#' See load_test_sensitivity() from
#' https://github.com/etchin/covid-testing/blob/master/R/testing.R for
#' code, original data, and correction to original data.
#'
#' @return a data frame with median, lower, and upper estimates by day
sens_by_day_data <- function() {
    structure(
        list(
            lower = c(
                0.38478627031159,
                0.57516731208943,
                0.604228529934369,
                0.615523437684205,
                0.634889417483735,
                0.658734818554505,
                0.683467990118063,
                0.704796690799587,
                0.715626316838783,
                0.708161673878986,
                0.677331542210053,
                0.628960600707915,
                0.57159750289502,
                0.513790902293819,
                0.46408945242676,
                0.427816321164893,
                0.397392733773658,
                0.362014429867096,
                0.318503776983481,
                0.294189654358031,
                0.268771607253158,
                0.243947269839712,
                0.219122932426265,
                0.194298595012819,
                0.169474257599372
            ),
            median = c(
                0.507191235810699,
                0.661649603931173,
                0.734764923759696,
                0.774741063381784,
                0.79627353568264,
                0.805680764747389,
                0.809281174661157,
                0.811418403300466,
                0.808536943707438,
                0.795106502715592,
                0.76729880640491,
                0.72809365784123,
                0.682172879336852,
                0.634218293204075,
                0.588911721755199,
                0.549417948006275,
                0.512833597788347,
                0.474738257636212,
                0.434289189502867,
                0.404954357014122,
                0.377078552659996,
                0.349468986554785,
                0.321859420449574,
                0.294249854344363,
                0.266640288239152
            ),
            upper = c(
                0.630862363616549,
                0.746928293064444,
                0.836882533155805,
                0.888289249824968,
                0.907958695940319,
                0.907512646103758,
                0.898572874917185,
                0.89014329958142,
                0.880756407692975,
                0.866326829447286,
                0.84380027533739,
                0.814246777046741,
                0.779767446556394,
                0.742463395847405,
                0.704435736900831,
                0.667378678469655,
                0.631358816394583,
                0.596035843288247,
                0.561457791759062,
                0.52922605439859,
                0.499840993181371,
                0.47051893746343,
                0.441196881745488,
                0.411874826027547,
                0.382552770309606
            )
        ),
        row.names = c(NA, -25L),
        class = "data.frame"
    )
}

gc2 <- function() {
    ## Quietly clear memory -- useful for GCP/AWS
    invisible(gc(FALSE))
}

## Summarizing helpers ----
summarize_sim_state <- function(sim_pop) {
    sim_pop %>%
        dplyr::group_by(time_step) %>%
        dplyr::summarize(
            n_susceptible = sum(State == 0),
            n_incubation_sub = sum(State == 1 & InfectionType == 0),
            n_infected_sub = sum(State == 2 & InfectionType == 0),
            n_incubation_clin = sum(State == 1 &
                                        InfectionType == 1),
            n_infected_clin = sum(State == 2 & InfectionType == 1),
            n_infected_all = sum(State %in% 1:2),
            n_infected_new = sum(NewInfection),
            n_recovered = sum(State == 3),
            n_quarantined = sum(DayInObs > 0),
            cume_n_detected = sum(DayOfDetection > 0),
            min_day_inf = min(DayOfInfection),
            max_day_inf = max(DayOfInfection),
            min_day_detect = min(DayOfDetection),
            max_day_detect = max(DayOfDetection)
        ) %>%
        dplyr::ungroup()
}

summarize_infectiousness <- function(sim_pop,
                                     if_weights,
                                     day_of_flight,
                                     subclin_infectious,
                                     min_if_weights = 0) {
    res <- tibble::tibble()
    for (w in min_if_weights) {
        x2 <- calculate_daily_infectiousness(
            sim_pop,
            if_weights,
            subclin_infectious,
            fast_mode = TRUE,
            min_if_weight = w
        )
        
        if (NROW(x2) > 0) {
            res <- dplyr::bind_rows(
                res,
                tibble::tibble(
                    time_step = unique(x2$time_step),
                    if_threshold = w,
                    n_infection_daily = sum(x2$n_infection_daily, na.rm = TRUE),
                    w_infection_daily = sum(x2$w_infection_daily, na.rm = TRUE),
                    n_active_infection = sum(x2$n_infection_daily > 0, na.rm = TRUE),
                    n_active_infection_subclin = sum(x2$n_infection_daily == .5, na.rm = TRUE)
                )
            )
        } else {
            ## If there are zero infections, calculate_infectiousness() returns
            ## an empty tibble().
            res <- dplyr::bind_rows(
                res,
                tibble::tibble(
                    time_step = unique(x2$time_step),
                    if_threshold = w,
                    n_infection_daily = 0,
                    w_infection_daily = 0,
                    n_active_infection = 0,
                    n_active_infection_subclin = 0
                )
            )
        }
        
    }
    res
}

## Analysis helpers ----
collect_results <- function(scenario_name) {
    all_files <- fs::dir_ls(dir_results(scenario_name),
                            recurse = TRUE,
                            glob = "*.RDS")
    
    purrr::map_df(.x = all_files,
                  .f = ~ readRDS(.x))
}

collect_and_munge_simulations <- function(scenario_name) {
    temp_x <- collect_results(scenario_name) %>%
        dplyr::filter(time_step >= 67,
                      !is.na(if_threshold)) %>%
        shift_time_steps(day_of_flight = 70) %>%
        add_cume_infections() %>%
        categorize_testing_types() %>%
        categorize_prob_inf() %>%
        categorize_sens_type() %>%
        categorize_prop_subclin() %>% 
        categorize_symptom_screening() %>%
        categorize_rt_multiplier() %>%
        categorize_risk_multiplier() %>%
        dplyr::mutate(sim_id = sprintf("%03d.%02d", round, rep)) %>%
        dplyr::select(
            testing_type,
            testing_cat,
            prob_inf,
            prob_inf_cat,
            sens_type,
            sens_cat,
            prop_subclin,
            prop_subclin_cat,
            symptom_screening,
            symptom_cat,
            risk_multiplier,
            risk_multi_cat,
            rapid_test_multiplier,
            rapid_test_cat,
            if_threshold,
            time_step,
            round,
            rep,
            sim_id,
            dplyr::everything()
        ) %>%
        dplyr::arrange(
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
            time_step
        )
    
    if (!tibble::has_name(temp_x, "n_test_false_pos") &
        !tibble::has_name(temp_x, "n_test_true_pos")) {
        temp_x <- temp_x %>%
            dplyr::mutate(n_test_false_pos = NA,
                          n_test_true_pos = NA)
    }
    
    temp_x
}

categorize_prop_subclin <- function(all_results) {
    all_results %>%
        mutate(prop_subclin_cat = factor(
            prop_subclin,
            levels = c(.3, .4),
            labels = c("30%", "40%"),
            ordered = TRUE
        ))
}

categorize_sens_type <- function(all_results) {
    all_results %>%
        mutate(sens_cat = factor(
            sens_type,
            levels = c("upper",
                       "median"),
            labels = c("Upper Bound",
                       "Median Value"),
            ordered = TRUE
        ))
}

categorize_testing_types <- function(all_results) {
    all_results %>%
        dplyr::mutate(testing_cat = factor(
            testing_type,
            levels = c(
                "no_testing",
                "no_testing_no_screening",
                "pcr_three_days_before",
                "pcr_three_days_before_5_day_quarantine_pcr",
                "rapid_test_same_day",
                "rapid_same_day_5_day_quarantine_pcr",
                "pcr_five_days_after",
                
                "pcr_two_days_before",
                "pcr_five_days_before",
                "pcr_seven_days_before",
                
                "pcr_three_days_before_7_day_quarantine_pcr",
                "pcr_three_days_before_14_day_quarantine_pcr",
                
                "pcr_two_days_before_5_day_quarantine_pcr",
                "pcr_five_days_before_5_day_quarantine_pcr",
                "pcr_seven_days_before_5_day_quarantine_pcr",
                
                "rapid_same_day_7_day_quarantine_pcr",
                "rapid_same_day_14_day_quarantine_pcr",
                
                "perfect_testing"
            ),
            labels = c(
                "No testing",
                "No testing, no screening",
                "PCR 3 days before",
                "PCR 3 days before + 5-day quarantine",
                "Same-day Rapid Test",
                "Same-day Rapid Test + 5-day quarantine",
                "PCR 5 days after",
                
                "PCR 2 days before",
                "PCR 5 days before",
                "PCR 7 days before",
                
                "PCR 3 days before + 7-day quarantine",
                "PCR 3 days before + 14-day quarantine",
                
                "PCR 2 days before + 5-day quarantine",
                "PCR 5 days before + 5-day quarantine",
                "PCR 7 days before + 5-day quarantine",
                
                "Same-day Rapid Test + 7-day quarantine",
                "Same-day Rapid Test + 14-day quarantine",
                
                "Daily perfect testing"
            ),
            ordered = TRUE
        ))
}

categorize_metric <- function(summarized_results) {
    summarized_results %>%
        dplyr::mutate(metric_cat =
                   factor(
                       metric,
                       levels = c(
                           "n_infection_daily",
                           "w_infection_daily",
                           "n_infected_all", 
                           "n_infected_new",
                           "n_active_infection", 
                           "n_active_infection_subclin", 
                           "n_active_infection_clin", 
                           "cume_n_infection_daily",
                           "cume_w_infection_daily",
                           "cume_n_infected_new",
                           "abs_cume_n_infected_new", 
                           "rel_n_infected_all", 
                           "abs_n_infection_daily",
                           "rel_n_infection_daily",
                           "ratio_n_infection_daily",
                           "abs_w_infection_daily",
                           "rel_w_infection_daily",
                           "ratio_w_infection_daily",
                           "abs_cume_n_infection_daily",
                           "rel_cume_n_infection_daily",
                           "ratio_cume_n_infection_daily",
                           "abs_cume_w_infection_daily",
                           "rel_cume_w_infection_daily",
                           "ratio_cume_w_infection_daily",
                           "frac_detected",
                           "frac_active_detected",
                           "any_positive_test",
                           "n_test_true_pos",
                           "n_test_false_pos",
                           "ratio_false_true",
                           "ratio_true_false",
                           "ppv",
                           "n_infected_day_of_flight",
                           "n_active_infected_day_of_flight",
                           "abs_n_active_infected_day_of_flight",
                           "rel_n_active_infected_day_of_flight",
                           "n_total_infections",
                           "abs_n_infected_all",
                           "abs_n_infected_new", 
                           "abs_n_active_infection",
                           "abs_n_active_infection_subclin", 
                           "abs_n_active_infection_clin",  
                           "rel_n_active_infection",
                           "rel_n_active_infection_subclin", 
                           "rel_n_active_infection_clin"
                       ),
                       labels = c(
                           "Infectious days",
                           "Weighted infections",
                           "All infections", 
                           "New infections",
                           "Active infections", 
                           "Active subclinical infections", 
                           "Active clinical infections",
                           "Cumulative infectious days",
                           "Cumulative weighted infections",
                           "Cumulative new infections",
                           "Abs. reduction in cumulative new infections", 
                           "Rel. reduction in total new infections", 
                           "Abs. difference in infectious days",
                           "Rel. difference in infectious days",
                           "Ratio of infectious days",
                           "Abs. difference in weighted infectiousness",
                           "Rel. difference in weighted infectiousness",
                           "Ratio of weighted infectiousness",
                           "Abs. difference in cumulative infection days",
                           "Rel. difference in cumulative infection days",
                           "Ratio of cumulative infection days",
                           "Abs. difference in cumulative infections",
                           "Rel. difference in cumulative infections",
                           "Ratio of cumulative infections",
                           "Fraction of infected detected",
                           "Fraction of active infections detected",
                           "Any positive test",
                           "True positives",
                           "False positives",
                           "False/true positive results",
                           "True/false positive results",
                           "Positive Predictive Value",
                           "Number infected on day of flight",
                           "Number active infections on day of flight",
                           "Abs. reduction in active infections on day of flight",
                           "Rel. reduction in active infections on day of flight",
                           "Total infections during observation",
                           "Abs. reduction in total infections",
                           "Abs. reduction in new infections", 
                           "Abs. reduction in active infections",
                           "Abs. reduction in active subclinical infections", 
                           "Abs. reduction in active clinical infections", 
                           "Rel. reduction in active infections",
                           "Rel. reduction in active subclinical infections", 
                           "Rel. reduction in active clinical infections"
                       ),
                       ordered = TRUE
                   ))
}

categorize_prob_inf <- function(all_results) {
    all_results %>%
        dplyr::mutate(prob_inf_cat =
                          factor(
                              prob_inf,
                              levels = c(50, 100, 200, 500, 1000, 
                                         1500, 2500, 5000) / 1000000,
                              labels = paste(c(5, 10, 20, 50, 
                                               100, 150, 250, 500),
                                             "per 100,000"),
                              ordered = TRUE
                          ))
}

categorize_symptom_screening <- function(all_results) {
    all_results %>%
        dplyr::mutate(symptom_cat =
                          factor(
                              symptom_screening,
                              levels = c(FALSE, TRUE),
                              labels = c("No symptom screening",
                                         "With symptom screening"),
                              ordered = TRUE
                          ))
}

categorize_rt_multiplier <- function(all_results) {
    all_results %>%
        dplyr::mutate(
            rapid_test_cat =
                factor(
                    rapid_test_multiplier,
                    levels = c(.6, .75, .9, 1), 
                    labels = sprintf("%i%%", round(c(.6, .75, .9, 1) * 100)),
                    ordered = TRUE
                )
        )
}

categorize_risk_multiplier <- function(all_results) {
    all_results %>%
        dplyr::mutate(risk_multi_cat =
                          factor(
                              risk_multiplier,
                              levels = c(1, 2, 4, 10),
                              labels = paste0(c(1, 2, 4, 10), "x"),
                              ordered = TRUE
                          ))
}

shift_time_steps <- function(all_results, day_of_flight = 70) {
    all_results %>%
        dplyr::mutate(relative_time = time_step - day_of_flight)
}

add_cume_infections <- function(all_results) {
    all_results %>%
        dplyr::arrange(
            testing_type,
            symptom_screening,
            risk_multiplier,
            rapid_test_multiplier,
            sens_type,
            prob_inf,
            if_threshold,
            round,
            rep,
            time_step
        ) %>%
        dplyr::group_by(
            testing_type,
            symptom_screening,
            risk_multiplier,
            rapid_test_multiplier,
            sens_type,
            prob_inf,
            round,
            rep,
            if_threshold
        ) %>%
        dplyr::mutate(
            cume_n_infection_daily = cumsum(n_infection_daily),
            cume_w_infection_daily = cumsum(w_infection_daily),
            cume_n_infected_new = cumsum(n_infected_new)
        ) %>%
        dplyr::ungroup()
}

summarize_results_column <- function(all_results, col) {
    all_results %>%
        dplyr::select({
            {
                col
            }
        }) %>%
        dplyr::mutate(
            n_total = dplyr::n(),
            n_infinite = sum(is.infinite({
                {
                    col
                }
            })),
            n_nan = sum(is.nan({
                {
                    col
                }
            })),
            n_missing = sum(is.na({
                {
                    col
                }
            }))
        ) %>%
        dplyr::filter(is.finite({
            {
                col
            }
        })) %>%
        dplyr::summarize(
            metric = rlang::as_label(dplyr::enquo(col)),
            n = dplyr::n(),
            mean = mean({
                {
                    col
                }
            }, na.rm = TRUE),
            sd = stats::sd({
                {
                    col
                }
            }, na.rm = TRUE),
            median = stats::median({
                {
                    col
                }
            }, na.rm = TRUE),
            p025 = stats::quantile({
                {
                    col
                }
            }, na.rm = TRUE, .025),
            p100 = stats::quantile({
                {
                    col
                }
            }, na.rm = TRUE, .10),
            p250 = stats::quantile({
                {
                    col
                }
            }, na.rm = TRUE, .25),
            p750 = stats::quantile({
                {
                    col
                }
            }, na.rm = TRUE, .75),
            p900 = stats::quantile({
                {
                    col
                }
            }, na.rm = TRUE, .90),
            p975 = stats::quantile({
                {
                    col
                }
            }, na.rm = TRUE, .975),
            min = min({
                {
                    col
                }
            }, na.rm = TRUE),
            max = max({
                {
                    col
                }
            }, na.rm = TRUE),
            n_missing = mean(n_missing),
            n_infinite = mean(n_infinite),
            n_nan = mean(n_nan),
            n_total = mean(n_total)
        ) %>%
        dplyr::ungroup()
}

## Triangle distribution code from Max Grossmann ----
## See: https://max.pm/posts/triangular_dist/
dtri <- function(x, a, b, c) {
    if (!all(a <= c && c <= b && a < b)) stop("It must be that a ≤ c ≤ b and a < b!");
    
    ifelse(x >= a & x <= b,
           ifelse(x <= c,
                  2*(x-a)/((b-a)*(c-a)),
                  2*(b-x)/((b-a)*(b-c))),
           0)
}

ptri <- function(x, a, b, c) {
    if (!all(a <= c && c <= b && a < b)) stop("It must be that a ≤ c ≤ b and a < b!");
    
    ifelse(x > a & x < b,
           ifelse(x <= c,
                  ((x-a)^2)/((b-a)*(c-a)),
                  1-((b-x)^2)/((b-a)*(b-c))),
           ifelse(x <= a,
                  0,
                  1))
}

qtri <- function(p, a, b, c) {
    if (!all(a <= c && c <= b && a < b)) stop("It must be that a ≤ c ≤ b and a < b!");
    
    ifelse(p > 0 & p < 1,
           ifelse(p <= ptri(c, a, b, c),
                  a+sqrt((a^2-a*b-(a-b)*c)*p),
                  b-sqrt(b^2-a*b+(a-b)*c+(a*b-b^2-(a-b)*c)*p)),
           NA)
}

rtri <- function(n, a, b, c) {
    if (!all(a <= c && c <= b && a < b)) stop("It must be that a ≤ c ≤ b and a < b!");
    
    qtri(stats::runif(n, min = 0, max = 1), a, b, c)
}
