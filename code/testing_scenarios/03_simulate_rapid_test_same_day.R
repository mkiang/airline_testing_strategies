## 02_simulate_rapid_test_same_day.R ----
##
## Simulate the same-day rapid antigen testing scenario,
##
## Each simulation is structured about the same. Use (the RStudio outline
## on the right (CMD+SHIFT+0) to jump around:
##      0. Set up
##      1. Create and burn in a simulated population (Days 0 - 62)
##      2. Pre-flight iterations (Days 63 to 69)
##      3. Day of flight iteration (Day 70).
##      4. Post-flight iterations (Days 71 to 84)
##
## In this scenario, during Step 3, we test the whole population with a
## rapid antigen test which has lower sensitivity.

simulate_rapid_test_same_day <- function(testing_type,
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
    ### Set file names ----
    sim_file <- return_sim_file_name(
        scenario_name = testing_type,
        symptom_screening = symptom_screening,
        prob_inf = prob_inf,
        prop_subclin = prop_subclin, 
        risk_multiplier = risk_multiplier,
        rapid_test_multiplier = rapid_test_multiplier,
        sens_type = sens_type, 
        round = round
    )
    
    ## Make sure there's a folder there
    fs::dir_create(dirname(sim_file))
    
    ### Start logging ----
    sink(sprintf("%s/log_%s.txt", dir_logs(), Sys.getpid()),
         append = TRUE)
    
    # Skip if we've done this
    if (!file_exists(sim_file)) {
        print_sim_start(
            testing_type,
            symptom_screening,
            prob_inf,
            sens_type,
            risk_multiplier,
            rapid_test_multiplier,
            prop_subclin, 
            round
        )
        
        ## Container to hold different results
        holder <- vector("list", n_reps)
        
        for (r in 1:n_reps) {
            ## 1. Create and burn-in simulation ----
            unpack_simulation_state(return_sim_state_file(round, r, prob_inf, prop_subclin))
            p_sens <- return_test_sensitivity(sens_type)
            
            print_sim_rep(r, days_incubation, days_symptomatic, prob_inf, prop_subclin)
            
            ## 2. Pre-flight iterations ----
            res <- tibble()
            for (k in n_burnin:(day_of_flight - 1)) {
                sim_pop <- increment_sim(
                    sim_pop = sim_pop,
                    t_step = k,
                    prob_inf = prob_inf,
                    if_weights = if_weights,
                    days_quarantine = days_quarantine
                )
                
                res <- bind_rows(
                    res,
                    left_join(
                        summarize_sim_state(sim_pop),
                        summarize_infectiousness(
                            sim_pop,
                            if_weights,
                            day_of_flight,
                            subclin_infectious
                        ),
                        by = "time_step"
                    )
                )
            }
            
            ## 3. Day of flight iteration ----
            ### Day of flight symptom screening (if applicable) ----
            if (symptom_screening) {
                sim_pop[DayInObs == 0 &
                            InfectionType == 1 &
                            State == 2, DayOfDetection := day_of_flight]
            }
            
            ### Perform same-day rapid antigen test ----
            sim_pop <- rapid_antigen_testing(
                sim_pop,
                t = day_of_flight,
                p_sens,
                p_spec,
                rapid_test_multiplier,
                days_incubation,
                days_symptomatic
            )
            
            ### Increment on the day of flight ----
            ## Keep this separate because day-of-flight has a risk multiplier
            sim_pop <- increment_sim(
                sim_pop = sim_pop,
                t_step = day_of_flight,
                prob_inf = prob_inf * risk_multiplier,
                if_weights = if_weights,
                days_quarantine = days_quarantine
            )
            
            res <- bind_rows(
                res,
                left_join(
                    summarize_sim_state(sim_pop) %>%
                        ## RECORD DAY OF FLIGHT TEST RESULTS
                        mutate(
                            n_test_false_pos = sum(sim_pop$sTested_fp),
                            n_test_true_pos = sum(sim_pop$sTested_tp)
                        ),
                    summarize_infectiousness(
                        sim_pop,
                        if_weights,
                        day_of_flight,
                        subclin_infectious
                    ),
                    by = "time_step"
                )
            )
            
            ## 4. Post-flight iterations ----
            for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
                ## We only want to count infections we could have potentially
                ## averted so after our t+3 (last possible testing day across
                ## all our scenarios), we stop infecting the population.
                if (k > (day_of_flight + 3)) {
                    sim_pop <- increment_sim(
                        sim_pop = sim_pop,
                        t_step = k,
                        prob_inf = 0,
                        if_weights = if_weights,
                        days_quarantine = days_quarantine
                    )
                } else {
                    sim_pop <- increment_sim(
                        sim_pop = sim_pop,
                        t_step = k,
                        prob_inf = prob_inf,
                        if_weights = if_weights,
                        days_quarantine = days_quarantine
                    )
                }
                
                res <- bind_rows(
                    res,
                    left_join(
                        summarize_sim_state(sim_pop),
                        summarize_infectiousness(
                            sim_pop,
                            if_weights,
                            day_of_flight,
                            subclin_infectious
                        ),
                        by = "time_step"
                        
                    )
                )
            }
            
            ## 5. Close out simulations ----
            ### Store summarized simulation results with parameter data ----
            rm(sim_pop)
            gc2()
            
            holder[[r]] <- res %>%
                mutate(
                    testing_type = testing_type,
                    symptom_screening = symptom_screening,
                    sens_type = sens_type,
                    round = round,
                    rep = r,
                    prob_inf = prob_inf,
                    prop_subclin = prop_subclin, 
                    risk_multiplier,
                    rapid_test_multiplier = rapid_test_multiplier
                )
        }
        
        ### Save summarized results ----
        saveRDS(bind_rows(holder),
                file = sim_file,
                compress = "xz")
        print(sprintf("Finished (%s): %s", Sys.time(), sim_file))
    } else {
        print(sprintf("Skipping (%s): %s", Sys.time(), sim_file))
    }
    
    ## Close log
    sink()
}
