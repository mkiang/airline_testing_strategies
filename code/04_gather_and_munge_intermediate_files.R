## 04_gather_and_munge_intermediate_files.R ----
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
source(here::here("code", "utils.R"))
fs::dir_create(here::here("data_raw"))

testing_scenarios <- basename(
    fs::dir_ls(here::here("intermediate_files"),
               type = "directory",
               regexp = "pcr|rapid|no_testing|perfect")
)

## Collect and save intermediate files in raw form ----
doParallel::registerDoParallel()
all_results <- foreach::foreach(s = testing_scenarios) %dopar% {
    s_name <- here::here("data_raw", sprintf("raw_simulations_%s.RDS", s))
    
    temp_x <- collect_and_munge_simulations(s)
    
    ### Save raw simulations ----
    saveRDS(temp_x,
            s_name,
            compress = "xz")
}
doParallel::stopImplicitCluster()
