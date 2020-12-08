simulate_perfect_daily_testing <- function(testing_type,
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
        testing_type,
        symptom_screening,
        prob_inf,
        sens_type,
        risk_multiplier,
        rapid_test_multiplier,
        round
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
            round
        )
        
        ## Container to hold different results
        holder <- vector("list", n_reps)
        
        for (r in 1:n_reps) {
            ## 1. Create and burn-in simulation ----
            unpack_simulation_state(return_sim_state_file(round, r, prob_inf))
            
            print_sim_rep(r, days_incubation, days_symptomatic, prob_inf)
            
            ## 2. Pre-flight iterations ----
            res <- tibble()
            for (k in n_burnin:(day_of_flight - 1)) {
                ## Test with perfect tests every day
                sim_pop <- pcr_test_subset_simpop(
                    sim_pop,
                    k,
                    n_subsets = 1,
                    subset_to_test = 1,
                    p_sens = return_test_sensitivity("perfect"),
                    p_spec = 1,
                    days_incubation,
                    days_symptomatic
                )
                
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
                        summarize_sim_state(sim_pop) %>%
                            mutate(
                                n_test_false_pos = sum(sim_pop$sTested_fp, na.rm = TRUE),
                                n_test_true_pos = sum(sim_pop$sTested_tp, na.rm = TRUE)
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
                
                sim_pop <- clear_old_tests(sim_pop)
            }
            
            ## 3. Day of flight iteration ----
            ### Day of flight symptom screening (if applicable)
            if (symptom_screening) {
                sim_pop[DayInObs == 0 &
                            InfectionType == 1 &
                            State == 2, DayOfDetection := day_of_flight]
            }
            
            ### Test with perfect tests every day ----
            sim_pop <- pcr_test_subset_simpop(
                sim_pop,
                k,
                n_subsets = 1,
                subset_to_test = 1,
                p_sens = return_test_sensitivity("perfect"),
                p_spec = 1,
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
                        mutate(
                            n_test_false_pos = sum(sim_pop$sTested_fp, na.rm = TRUE),
                            n_test_true_pos = sum(sim_pop$sTested_tp, na.rm = TRUE)
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
            
            sim_pop <- clear_old_tests(sim_pop)
            
            ## 4. Post-flight iterations ----
            for (k in (day_of_flight + 1):(day_of_flight + n_outcome_day)) {
                ## Test with perfect tests every day
                sim_pop <- pcr_test_subset_simpop(
                    sim_pop,
                    k,
                    n_subsets = 1,
                    subset_to_test = 1,
                    p_sens = return_test_sensitivity("perfect"),
                    p_spec = 1,
                    days_incubation,
                    days_symptomatic
                )
                
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
                        summarize_sim_state(sim_pop) %>%
                            mutate(
                                n_test_false_pos = sum(sim_pop$sTested_fp, na.rm = TRUE),
                                n_test_true_pos = sum(sim_pop$sTested_tp, na.rm = TRUE)
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
                
                sim_pop <- clear_old_tests(sim_pop)
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

## Close connections
doParallel::stopImplicitCluster()
