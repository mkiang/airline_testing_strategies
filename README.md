
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Routine testing strategies for airline travel during the COVID-19 pandemic: a simulation analysis

<p align="center">
<img src="./plots/fig1_cume_inf_over_time.jpg" width="600px" style="display: block; margin: auto;" />
</p>

## Introduction

Reproducible code for our paper, [*Routine testing strategies for
airline travel during the COVID-19 pandemic: a simulation
analysis*](TODO), in which we use microsimulation to estimate the
effectiveness of different testing strategies on reducing SARS-CoV-2
transmission in a cohort of hypothetical airline travels.

The full citation is:

> Kiang MV, Chin ET, Huynh BQ, Chapman LAC, Rodríguez-Barraquer I,
> Greenhouse B, Rutherford GW, Bibbins-Domingo K, Havlir D, Basu S, and
> Lo NC. Routine asymptomatic testing strategies for airline travel
> during the COVID-19 pandemic: a simulation analysis. 2020.

### Abstract

**Background:** Airline travel has been significantly reduced during the
COVID-19 pandemic due to concern for individual risk of SARS-CoV-2
infection and population-level transmission risk from importation.
Routine viral testing strategies for COVID-19 may facilitate safe
airline travel through reduction of individual and/or population-level
risk, although the effectiveness and optimal design of these
“test-and-travel” strategies remain unclear.

**Methods:** We developed a microsimulation of SARS-CoV-2 transmission
in a cohort of airline travelers to evaluate the effectiveness of
various testing strategies to reduce individual risk of infection and
population-level risk of transmission. We evaluated five testing
strategies in asymptomatic passengers: i) anterior nasal polymerase
chain reaction (PCR) within 3 days of departure; ii) PCR within 3 days
of departure and PCR 5 days after arrival; iii) rapid antigen test on
the day of travel (assuming 90% of the sensitivity of PCR during active
infection); iv) rapid antigen test on the day of travel and PCR 5 days
after arrival; and v) PCR within 3 days of arrival alone. The travel
period was defined as three days prior to the day of travel and two
weeks following the day of travel, and we assumed passengers followed
guidance on mask wearing during this period. The primary study outcome
was cumulative number of infectious days in the cohort over the travel
period (population-level transmission risk); the secondary outcome was
the proportion of infectious persons detected on the day of travel
(individual-level risk of infection). Sensitivity analyses were
conducted.

**Findings:** Assuming a community SARS-CoV-2 incidence of 50 daily
infections, we estimated that in a cohort of 100,000 airline travelers
followed over the travel period, there would be a total of 2,796 (95%
UI: 2,031, 4,336) infectious days with 229 (95% UI: 170, 336) actively
infectious passengers on the day of travel. The pre-travel PCR test
(within 3 days prior to departure) reduced the number of infectious days
by 35% (95% UI: 27, 42) and identified 88% (95% UI: 76, 94) of the
actively infectious travelers on the day of flight; the addition of PCR
5 days after arrival reduced the number of infectious days by 79% (95%
UI: 71, 84). The rapid antigen test on the day of travel reduced the
number of infectious days by 32% (95% UI: 25, 39) and identified 87%
(95% UI: 81, 92) of the actively infectious travelers; the addition of
PCR 5 days after arrival reduced the number of infectious days by 70%
(95% UI: 65, 75). The post-travel PCR test alone (within 3 days of
landing) reduced the number of infectious days by 42% (95% UI: 31, 51).
The ratio of true positives to false positives varied with the incidence
of infection. The overall study conclusions were robust in sensitivity
analysis.

**Interpretation:** Routine asymptomatic testing for COVID-19 prior to
travel can be an effective strategy to reduce individual risk of
COVID-19 infection during travel, although post-travel testing with
abbreviated quarantine is likely needed to reduce population-level
transmission due to importation of infection when traveling from a high
to low incidence setting.

### Issues

Please report issues via email or the [issues
page](https://github.com/mkiang/airline_testing_strategies/issues).

## Structure

This project is structured as follows:

-   `code`: Contains all code used for this project. Designed to be run
    sequentially. A brief overview of the code files is provided below.
-   `data`: Contains all data necessary to reproduce our tables and
    plots.
-   `intermediate_files`: Contains temporary files not necessary to
    reproduce our tables and plots.
-   `output`: Contains numerical representations (in `csv` form) of all
    figures.
-   `plots`: Contains all figures used the manuscript in `jpg` and
    `pdf`.
-   `rmds`: Contains useful rmarkdown files such as reproducing our
    tables and providing session information.

Note, we use the `config.yml` file to modify high-level project-wide
parameters.

### Reproducibility

All information necessary for full reproducibility, including package
version numbers, is available in `./rmds/02_session_info.html`. This
project uses the
[`renv`](https://rstudio.github.io/renv/articles/renv.html) package for
package version control. To use this, open the project in
[RStudio](https://rstudio.com/) and run `renv::restore()`.

### `config.yml`

Simulation parameters such as the length of quarantine, number of
burn-in days, number of passengers, proportion of subclinical
infections, daily probability of infection, etc. can be changed in the
`config.yml` file. In addition, the `config.yml` file lets you specify
different sensitivity analysis configurations and allows you to vary one
(or more) parameters across the default parameter sweep. See the
`config.yml` file for details.

## `./code`

The analytic pipeline in the `./code`/ folder is designed to have each
file run sequentially. That is, each file performs a discrete task and
some tasks may be dependent on the output of previous tasks.

-   `utils.R`: Contains almost all helper functions and functions
    necessary to run the simulations or create the plots. Most functions
    within `utils.R` uses the
    [`roxygen2`](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html)
    documentation template. See this documentation for expected input,
    output, and notes about each function.
-   `01_create_simulated_populations.R`: Creates and performs the
    burn-in for all simulated passenger cohorts. In addition, each
    simulation state is saved with the random seed used to create the
    simulation as well as the seed necessary to resume the simulation.
    Each testing scenario calls one of the simulation states created
    here (using `unpack_simulation_state()`).
-   `02_run_main_simulations.R`: Run all testing scenarios used for the
    main results. Each simulation result is saved in
    `./intermediate_files` to be collected later using
    `04_summarize_infection_quantities.R` and
    `05_summarize_testing_quantites.R` and then saved in `./data`.
-   `03_run_sensitivity_analyses.R`: Run all testing scenarios used for
    the main results. Each simulation result is saved in
    `intermediate_files` to be collected later using
    `04_summarize_infection_quantities.R` and
    `05_summarize_testing_quantites.R` and then saved in `./data`.
-   `04_summarize_infection_quantities.R`: Collects all the files
    created in `02` and `03` in `./intermediate_files`, processes them,
    and then saves the result in `./data`. This file saves the
    time-varying quantities of interest and our primary end points.
-   `05_summarize_testing_quantites.R`: Collects all the files created
    in `02` and `03` in `./intermediate_files`, processes them, and then
    saves the result in `./data`. This file saves the time-invariant
    quantities (e.g., testing results).

Other files can be run in any order after `04` and `05` have
successfully run. These files create plots used in the main manuscript
or associated appendix.

### `./code/testing_scenarios`

Testing scenarios are written as functions that take in all necessary
parameters and return no output, instead they save their intermediate
results in a temporary location. The structure each function is the
same:

1.  Load a simulated population using `unpack_simulation_state()`
2.  Increment pre-flight time steps
3.  Increment the day-of-flight time step
4.  Increment post-flight time steps
5.  Repeat these steps `n_reps` number of times
6.  Save all `n_reps` realizations into a temporary file in
    `./intermediate_files`

The scenarios themselves change in what happens within each stage above.
For example, pre-flight testing strategies would occur in Step 2. The
day-of-flight risk multiplier only impacts Step 3. Post-flight testing
strategies occur in Step 4. See the documentation within each scenario
file in `./code/testing_scenarios` for details.

## Authors (alphabetical):

-   [Sanjay
    Basu](https://primarycare.hms.harvard.edu/faculty-staff/sanjay-basu)
    (![Github](http://i.imgur.com/9I6NRUm.png):
    [sanjaybasu](https://github.com/sanjaybasu) \|
    ![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@sanjaybmdphd](https://twitter.com/sanjaybmdphd))
-   [Kirsten
    Bibbins-Domingo](https://profiles.ucsf.edu/kirsten.bibbins-domingo)
    (![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@KBibbinsDomingo](https://twitter.com/KBibbinsDomingo))
-   [Lloyd Chapman](https://profiles.ucsf.edu/lloyd.chapman)
    (![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@lloyd\_chapman\_](https://twitter.com/lloyd_chapman_))
-   Elizabeth Chin (![Github](http://i.imgur.com/9I6NRUm.png):
    [etchin](https://github.com/etchin) \|
    ![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@elizabethtchin](https://twitter.com/elizabethtchin))
-   [Bryan Greenhouse]()
-   [Diane Havlir]()
-   [Benjamin Huynh](https://benhuynh.github.io/)
    (![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@benqhuynh](https://twitter.com/benqhuynh))
-   [Mathew Kiang](https://mathewkiang.com)
    (![Github](http://i.imgur.com/9I6NRUm.png):
    [mkiang](https://github.com/mkiang) \|
    ![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@mathewkiang](https://twitter.com/mathewkiang))
-   [Nathan Lo](https://profiles.ucsf.edu/nathan.lo)
    (![Github](http://i.imgur.com/9I6NRUm.png):
    [NathanLo3](https://github.com/NathanLo3) \|
    ![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@NathanLo3579](https://twitter.com/NathanLo3579))
-   [Isabel
    Rodríguez-Barraquer](https://profiles.ucsf.edu/isabel.rodriguez-barraquer)
    (![Github](http://i.imgur.com/9I6NRUm.png):
    [isabelrodbar](https://github.com/isabelrodbar) \|
    ![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@isabelrodbar](https://twitter.com/isabelrodbar))
-   [George Rutherford](https://profiles.ucsf.edu/george.rutherford)
    (![Twitter](http://i.imgur.com/wWzX9uB.png):
    [@rutherford\_ucsf](https://twitter.com/rutherford_ucsf))

## Attribution

This simulation code is heavily based on the paper, [*Frequency of
Routine Testing for Coronavirus Disease 2019 (COVID-19) in High-risk
Healthcare Environments to Reduce
Outbreaks*](https://academic.oup.com/cid/advance-article/doi/10.1093/cid/ciaa1383/5939986),
by Chin et al. and the accompanying [Github
repository](https://github.com/etchin/covid-testing).
