# Script for running scenarios for the RIVM Scenarios Project
# This script runs simulations for Wave 1

# preamble ---------------------------------------------------------
# This script will load necessary packages, data sources, run 
# scenarios.
# ------------------------------------------------------------------

# Options ----------------------------------------------------------
# suppress dplyr::summarise() warnings
options(dplyr.summarise.inform = FALSE)

# Load required packages -------------------------------------------
library(deSolve)
#library(reshape2)
library(tidyverse)
library(readxl)
#library(rARPACK)
library(lubridate)
library(foreach)
library(doParallel)
library(here)

# Source required files --------------------------------------------
source("R/na_to_zero.R")
source("R/age_struct_seir_simple.R")
source("R/postprocess_age_struct_model_output_simple.R")
source("R/summarise_results_simple.R")
source("R/get_foi_simple.R")
source("R/choose_contact_matrix_simple.R")
# -------------------------------------------------------------------
# Define population size --------------------------------------------
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 
              0.13083463,0.14514332, 0.12092904, 0.08807406, 
              0.04622194)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

# ve estimates ------------------------------------------------------
ve <- read_excel("inst/extdata/inputs/ve_estimates/ve_dat.xlsx", sheet = "wildtype") %>%
  group_by(dose, age_group, outcome) %>%
  summarise(mean_ve = mean(ve))

ve_inf <- ve %>% 
  filter(outcome == "infection",
         dose == "d2")

ve_hosp <- ve %>% 
  filter(outcome == "hospitalisation",
         dose == "d2")

ve_trans <- ve %>% 
  filter(outcome == "transmission",
         dose == "d2")

# probabilities -------------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/inputs/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_hospital2death <- c(rep(0, 5), 0.01, 0.04, 0.12, 0.29) # (after ICU)
# p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien

# delays --------------------------------------------------------------
time_symptom2admission <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC <- 2.28
time_IC2hospital <- 15.6
time_hospital2discharge <- 10.1 # (after ICU)
time_admission2death <- 7
time_IC2death <- 19
time_hospital2death <- 10 # (after ICU)

# define transition rates ---------------------------------------------
time_in_i <- ((1-p_infection2admission) * 2) + (p_infection2admission * time_symptom2admission)
i2r    <- (1-p_infection2admission) / time_in_i          
i2h    <- p_infection2admission / time_in_i 

time_in_iv <- ((1-(ve_hosp$mean_ve*p_infection2admission)) * 2) + ((ve_hosp$mean_ve*p_infection2admission) * time_symptom2admission)
iv2rv    <- (1-(ve_hosp$mean_ve*p_infection2admission)) / time_in_iv          
iv2hv    <- (ve_hosp$mean_ve*p_infection2admission) / time_in_iv 

time_in_h <- (p_admission2IC * time_admission2IC) + (p_admission2death * time_admission2death) + ((1 - (p_admission2IC + p_admission2death)) * time_admission2discharge)
h2ic   <- p_admission2IC / time_in_h                    
h2d    <- p_admission2death / time_in_h  
h2r    <- (1 - (p_admission2IC + p_admission2death)) / time_in_h

time_in_ic <- (p_IC2hospital * time_IC2hospital) + ((1 - p_IC2hospital) * time_IC2death)
ic2hic <- p_IC2hospital / time_in_ic
ic2d   <- (1 - p_IC2hospital) / time_in_ic

time_in_hic <- (p_hospital2death * time_hospital2death) + ((1 - p_hospital2death) * time_hospital2discharge)
hic2d  <- p_hospital2death / time_in_hic         
hic2r  <- (1 - p_hospital2death) / time_in_hic

# contact matrices --------------------------------------------------
path <- "/rivm/s/ainsliek/data/contact_matrices/converted/"
# path <- "inst/extdata/inputs/contact_matrices/converted/"
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
june_2022      <- readRDS(paste0(path,"transmission_matrix_june_2022.rds"))

# specify initial model parameters ---------------------------------
# parameters must be in a named list
params <- list(N = n_vec,  # population size
               # rates
               beta = 0.0005,
               beta1 = 0, #0.14,
               sigma = 0.5,
               epsilon = 0.00,
               omega = 0.0038, # 60% immunity after 8 months from Exp dist
               gamma = i2r,
               gamma_v = iv2rv,
               h = i2h,
               h_v = iv2hv,
               i1 = h2ic,
               d = h2d,
               r = h2r,
               i2 = ic2hic,
               d_ic = ic2d,
               d_hic = hic2d,
               r_ic = hic2r,
               # simulation start time
               calendar_start_date = as.Date("2020-01-01"), 
               # contact matrices for different levels of NPIs
               c_start = june_2022$mean, # use Pico 8 (when available)
               c_lockdown = april_2020$mean,
               c_open = april_2020$mean, # keep lockdown contact matrix
               keep_cm_fixed = FALSE,
               # IC admission thresholds
               thresh_o = 1,
               thresh_l = 40,
               # vaccination parameters
               alpha = c(rep(0,9)),
               t_vac_start = NULL,
               t_vac_end = NULL,
               eta = 1, 
               eta_trans = 1
)

# Specify initial conditions --------------------------------------
# get initial conditions from last time point (275) of Wave 1
# read in results from Wave 1
scenarioA_wave1 <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave1_scenarioA.rds")
scenarioB_wave1 <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave1_scenarioB.rds")
scenarioC_wave1 <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave1_scenarioC.rds")
scenarioD_wave1 <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave1_scenarioD.rds")
scenarioE_wave1 <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave1_scenarioE.rds")

# get tail of dataset for each simulation
initA <- lapply(scenarioA_wave1, tail, 1)
initB <- lapply(scenarioB_wave1, tail, 1)
initC <- lapply(scenarioC_wave1, tail, 1)
initD <- lapply(scenarioD_wave1, tail, 1)
initE <- lapply(scenarioE_wave1, tail, 1)

# Run forward simulations --------------------------------------------
# Scenario A: no measures
# Scenario B: voluntary
# Scenario C: R < 1 @ high inf rate
# Scenario D: R < 1 @ low inf rate
# Scenario E: zero COVID
t_start <- 0 #as.numeric(initA[[1]][1])
t_end <- t_start + 270
times <- as.integer(seq(t_start, t_end, by = 1))
# betas <- readRDS("../vacamole/inst/extdata/results/model_fits/beta_draws.rds")
# # sample 100 betas from last time window
# betas100 <- sample(betas[[1]]$beta, 100)

# register parallel backend
registerDoParallel(cores=15)
n_sim <- 100
# Scenario A: no measures ----
scenarioA <- foreach(i = 1:n_sim) %dopar% {
  flag_open <- 0
  init <- unlist(initA[[i]][-1])
  init[1] <- t_start
  # params$beta <- betas100[i] * 1.5
  params$c_start <- june_2022[[i]]
  params$keep_cm_fixed <- TRUE # force contact matrix to stay fixed
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, times, age_struct_seir_simple, params, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioA, "/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioA.rds")
doParallel::stopImplicitCluster()

# Scenario B: voluntary ----
# registerDoParallel(cores=15)
scenarioB <- foreach(i = 1:n_sim) %dopar% {
  flag_open <- 0
  init <- unlist(initB[[i]][-1])
  init[1] <- t_start
  params$keep_cm_fixed <- FALSE
  #params$beta <- betas100[i] * 1.5
  params$c_start <- june_2022[[i]]
  params$c_lockdown <- june_2020[[i]]
  params$c_open <- params$c_lockdown
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, times, age_struct_seir_simple, params, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioB, "/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioB.rds")
doParallel::stopImplicitCluster()

# Scenario C: R<1 @ LOW incidence ----
registerDoParallel(cores=15)
scenarioC <- foreach(i = 1:n_sim) %dopar% {
  flag_open <- 0
  init <- unlist(initC[[i]][-1])
  init[1] <- t_start
  params$keep_cm_fixed <- FALSE
  # params$beta <- betas100[i] * 1.5
  params$c_start <- june_2022[[i]]
  params$c_lockdown <- april_2020[[i]]
  params$c_open <- params$c_lockdown
  params$thresh_l <- 10
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, times, age_struct_seir_simple, params, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioC, "/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioC.rds")
doParallel::stopImplicitCluster()

# Scenario D: R<1 @ HIGH incidence ----
registerDoParallel(cores=15)
scenarioD <- foreach(i = 1:n_sim) %dopar% {
  flag_open <- 0
  init <- unlist(initD[[i]][-1])
  init[1] <- t_start
  params$keep_cm_fixed <- FALSE
  # params$beta <- betas100[i] * 1.5
  params$c_start <- june_2022[[i]]
  params$c_lockdown <- april_2020[[i]]
  params$c_open <- params$c_lockdown
  params$thresh_l <- 40
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, times, age_struct_seir_simple, params, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioD, "/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioD.rds")
doParallel::stopImplicitCluster()

# Scenario E: zero covid ----
registerDoParallel(cores=15)
scenarioE <- foreach(i = 1:n_sim) %dopar% {
  flag_open <- 0
  init <- unlist(initE[[i]][-1])
  init[1] <- t_start
  params$keep_cm_fixed <- FALSE
  # params$beta <- betas100[i] * 1.5
  params$c_start <- june_2022[[i]]
  params$c_lockdown <- april_2020[[i]]
  params$c_open <- params$c_lockdown
  params$thresh_l <- 1
  params$thresh_o <- 0
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, times, age_struct_seir_simple, params, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioE, "/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioE.rds")
doParallel::stopImplicitCluster()

# Post-process scenario runs ---------------------------------------------------
# Results must be in a csv file that contains only the following columns (in any
# order). No additional columns are allowed.
# - origin_date (date):	Date as YYYY-MM-DD, last day (Monday) of submission window
# - scenario_id	(string):	A specified "scenario ID"
# - target_variable	(string):	"inc case", "inc death", "inc hosp", "inc icu", "inc infection"
# - horizon	(string):	The time horizon for the projection relative to the origin_date (e.g., X wk)
# - target_end_date	(date):	Date as YYYY-MM-DD, the last day (Saturday) of the target week
# - location (string): An ISO-2 country code or "H0" ("NL" for the Netherlands)
# - sample	(numeric):	A integer corresponding to the sample index (used to plot trajectories)
# - value	(numeric):	The projected count, a non-negative integer number of new cases or deaths in the epidemiological week

# wrangle Scenario A output ----------------------------------------------------
# read in saved output from model runs
scenarioA <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioA.rds")
scenarioB <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioB.rds")
scenarioC <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioC.rds")
scenarioD <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioD.rds")
scenarioE <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave2_scenarioE.rds")

# create empty lists to store outputs
outA <- list()
outB <- list()
outC <- list()
outD <- list()
outE <- list()

# loop over samples and summarize results for each scenario
for(s in 1:n_sim){
  # specify shared parameter values (transmission rate and starting contact matrix)
  #params$beta <- betas100[s] * 1.5
  params$c_start <- june_2022[[s]] # Pico 8 when available
  
  # Scenario A - no measures
  params$keep_cm_fixed <- TRUE
  seir_outputA <- postprocess_age_struct_model_output_simple(scenarioA[[s]])
  seir_outcomesA <- summarise_results_simple(seir_outputA, params = params, t_vec = times) %>%
    mutate(sample = s)
  outA[[s]] <- seir_outcomesA
  
  # Scenario B - voluntary
  params$keep_cm_fixed <- FALSE
  params$c_lockdown <- june_2020[[s]]
  params$c_open <- params$c_lockdown
  
  seir_outputB <- postprocess_age_struct_model_output_simple(scenarioB[[s]])
  seir_outcomesB <- summarise_results_simple(seir_outputB, params = params, t_vec = times) %>%
    mutate(sample = s)
  outB[[s]] <- seir_outcomesB
  
  # Scenario C - R<1 @ low incidence
  params$keep_cm_fixed <- FALSE
  params$c_lockdown <- april_2020[[s]]
  params$c_open <- params$c_lockdown
  params$thresh_l <- 10
  
  seir_outputC <- postprocess_age_struct_model_output_simple(scenarioC[[s]])
  seir_outcomesC <- summarise_results_simple(seir_outputC, params = params, t_vec = times) %>%
    mutate(sample = s)
  outC[[s]] <- seir_outcomesC
  
  # Scenario D - R<1 @ high incidence
  params$keep_cm_fixed <- FALSE
  params$c_lockdown <- april_2020[[s]]
  params$c_open <- params$c_lockdown
  params$thresh_l <- 40
  
  seir_outputD <- postprocess_age_struct_model_output_simple(scenarioD[[s]])
  seir_outcomesD <- summarise_results_simple(seir_outputD, params = params, t_vec = times) %>%
    mutate(sample = s)
  outD[[s]] <- seir_outcomesD
  
  # Scenario E - zero covid
  params$keep_cm_fixed <- FALSE
  params$c_lockdown <- april_2020[[s]]
  params$c_open <- params$c_lockdown
  params$thresh_l <- 1
  params$thresh_o <- 0
  
  seir_outputE <- postprocess_age_struct_model_output_simple(scenarioE[[s]])
  seir_outcomesE <- summarise_results_simple(seir_outputE, params = params, t_vec = times) %>%
    mutate(sample = s)
  outE[[s]] <- seir_outcomesE
  
}

# create data frames
dfA <- bind_rows(outA) %>% mutate(scenario_id = "A") 
dfB <- bind_rows(outB) %>% mutate(scenario_id = "B") 
dfC <- bind_rows(outC) %>% mutate(scenario_id = "C")
dfD <- bind_rows(outD) %>% mutate(scenario_id = "D")
dfE <- bind_rows(outE) %>% mutate(scenario_id = "E")

# create single post-processed data frame to save to directory
df_wave2 <- bind_rows(dfA, dfB, dfC, dfD, dfE)

saveRDS(df_wave2, "inst/extdata/results/wave2_results.rds")

# free unused memory
gc()








