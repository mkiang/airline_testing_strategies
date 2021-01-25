## Testing scenarios ----
## These functions are held as their own files. See each file for details. 
## Null models
source(here::here("code", "testing_scenarios", "00a_simulate_null_model.R"))
source(here::here("code", "testing_scenarios", "00b_simulate_perfect_daily_testing.R"))
## Main testing scenarios
source(here::here("code", "testing_scenarios", "01_simulate_pcr_three_days_before.R"))
source(here::here("code", "testing_scenarios", "02_simulate_pcr_3_days_before_5_day_quarantine_pcr.R"))
source(here::here("code", "testing_scenarios", "03_simulate_rapid_test_same_day.R"))
source(here::here("code", "testing_scenarios", "04_simulate_rapid_antigen_same_day_5_day_quarantine_pcr.R"))
source(here::here("code", "testing_scenarios", "05_simulate_pcr_five_days_after.R"))
## Sensitivity testing scenarios
source(here::here("code", "testing_scenarios", "01a_simulate_pcr_two_days_before.R"))
source(here::here("code", "testing_scenarios", "01b_simulate_pcr_five_days_before.R"))
source(here::here("code", "testing_scenarios", "01c_simulate_pcr_seven_days_before.R"))
source(here::here("code", "testing_scenarios", "02a_simulate_pcr_3_days_before_7_day_quarantine_pcr.R"))
source(here::here("code", "testing_scenarios", "02b_simulate_pcr_3_days_before_14_day_quarantine_pcr.R"))
source(here::here("code", "testing_scenarios", "02c_simulate_pcr_2_days_before_5_day_quarantine_pcr.R"))
source(here::here("code", "testing_scenarios", "02d_simulate_pcr_5_days_before_5_day_quarantine_pcr.R"))
source(here::here("code", "testing_scenarios", "02e_simulate_pcr_7_days_before_5_day_quarantine_pcr.R"))
source(here::here("code", "testing_scenarios", "04a_simulate_rapid_antigen_same_day_7_day_quarantine_pcr.R"))
source(here::here("code", "testing_scenarios", "04b_simulate_rapid_antigen_same_day_14_day_quarantine_pcr.R"))

#' Wrapper function for selecting testing scenarios
#' 
#' Note, this is just a wrapper for passing things along to the underlying
#' scenario functions. See each scenario function for details. 
#'
#' @param testing_type testing scenario short code
#' @param prob_inf probability of infection 
#' @param sens_type sensitivity type for return_test_sensitivity()
#' @param risk_multiplier day of travel risk multiplier
#' @param rapid_test_multiplier discount for rapid test sensitivity relative to PCR
#' @param symptom_screening TRUE/FALSE for day of flight symptomatic screening
#' @param round which round of simulation is this?
#' @param n_reps number of repetitions per round of simulation
#' @param prop_subclin proportion of infections that are subclinical
#' @param subclin_infectious discount for subclinical infectiousness
#' @param n_burnin number of days to burn in
#' @param day_of_flight day of travel 
#' @param n_outcome_day days post-travel to follow the cohort
#' @param days_quarantine duration of quarantine
#'
#' @return none (saves file into intermediate_files folder)
run_and_save_simulation <- function(testing_type,
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
                                    days_quarantine) {
    if (!(testing_type %in% c("no_testing",
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
    )
    )) {
        stop("Not a valid testing_type")
    }
    
    if (testing_type == "no_testing") {
        simulate_null_model(
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
    
    if (testing_type == "perfect_testing") {
        simulate_perfect_daily_testing(
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
    
    if (testing_type == "pcr_three_days_before") {
        simulate_pcr_three_days_before(
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
    
    if (testing_type == "pcr_three_days_before_5_day_quarantine_pcr") {
        simulate_pcr_before_5_day_quarantine_pcr(
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
    
    if (testing_type == "rapid_test_same_day") {
        simulate_rapid_test_same_day(
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
    
    if (testing_type == "rapid_same_day_5_day_quarantine_pcr") {
        simulate_rapid_antigen_same_day_5_day_quarantine_pcr(
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
    
    if (testing_type == "pcr_five_days_after") {
        simulate_pcr_five_days_after(
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
    
    if (testing_type == "pcr_two_days_before") {
        simulate_pcr_two_days_before(
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
    
    if (testing_type == "pcr_five_days_before") {
        simulate_pcr_five_days_before(
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
    
    if (testing_type == "pcr_seven_days_before") {
        simulate_pcr_seven_days_before(
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
    
    if (testing_type == "pcr_three_days_before_7_day_quarantine_pcr") {
        simulate_pcr_before_7_day_quarantine_pcr(
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
    
    if (testing_type == "pcr_three_days_before_14_day_quarantine_pcr") {
        simulate_pcr_before_14_day_quarantine_pcr(
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
    
    if (testing_type == "pcr_two_days_before_5_day_quarantine_pcr") {
        simulate_pcr_two_before_5_day_quarantine_pcr(
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
    
    if (testing_type == "pcr_five_days_before_5_day_quarantine_pcr") {
        simulate_pcr_five_before_5_day_quarantine_pcr(
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
    
    if (testing_type == "pcr_seven_days_before_5_day_quarantine_pcr") {
        simulate_pcr_seven_before_5_day_quarantine_pcr(
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
    
    if (testing_type == "rapid_same_day_7_day_quarantine_pcr") {
        simulate_rapid_antigen_same_day_7_day_quarantine_pcr(
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
    
    if (testing_type == "rapid_same_day_14_day_quarantine_pcr") {
        simulate_rapid_antigen_same_day_14_day_quarantine_pcr(
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
}
