# Script for running scenarios for the Eupropean Scenario Hub
# URL: https://github.com/covid19-forecast-hub-europe/covid19-scenario-hub-europe#readme

# preamble ---------------------------------------------------------
# This script will load necessary packages, data sources, fit the model to data,
# and then run scenarios.
# All scenarios will be using Dutch data
# TODO: generalize to other European countries
# TODO: break up code into work chunks
# ------------------------------------------------------------------

# Options ----------------------------------------------------------
# suppress dplyr::summarise() warnings
options(dplyr.summarise.inform = FALSE)

# Load required packages/functions ---------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(rARPACK)
library(readr)
library(lubridate)
library(foreach)
library(doParallel)
library(here)

# source("R/convert_vac_schedule2.R")
source("R/na_to_zero.R")
# source("R/calc_waning.R")
source("R/age_struct_seir_simple.R")
source("R/postprocess_age_struct_model_output2.R")
source("R/summarise_results.R")
source("R/get_foi.R")
# -------------------------------------------------------------------
# Define population size --------------------------------------------
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 
              0.13083463,0.14514332, 0.12092904, 0.08807406, 
              0.04622194)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

# probabilities -------------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/inputs/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_hospital2death <- c(rep(0, 5), 0.01, 0.04, 0.12, 0.29) # (after ICU)
p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 
                       0.409) # from Jantien

# delays --------------------------------------------------------------
time_symptom2admission <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 
                            4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC <- 2.28
time_IC2hospital <- 15.6
time_hospital2discharge <- 10.1 # (after ICU)
time_admission2death <- 7
time_IC2death <- 19
time_hospital2death <- 10 # (after ICU)

# define transition rates ---------------------------------------------
i2r    <- (1-p_infection2admission) / 2                    # I -> R
i2h    <- p_infection2admission / time_symptom2admission   # I -> H

h2ic   <- p_admission2IC / time_admission2IC               # H -> IC
h2d    <- p_admission2death / time_admission2death         # H -> D
h2r    <- (1 - (p_admission2IC + p_admission2death)) / time_admission2discharge
# H -> R

ic2hic <- p_IC2hospital / time_IC2hospital                 # IC -> H_IC
ic2d   <- (1 - p_IC2hospital) / time_IC2death              # IC -> D

hic2d  <- p_hospital2death / time_hospital2death           # H_IC -> D
hic2r  <- (1 - p_hospital2death) / time_hospital2discharge # H_IC -> R

# determine waning rate from Erlang distribution --------------------
# We want the rate that corresponds to a 60% reduction in immunity after 
#   - 3 months (92 days) or
#   - 8 months (244 days)

# we need to solve the following equation for lambda (waning rate)
# tau = time since recovery
# p = probability still immune
Fk <- function(lambda, tau, p){
  exp(-tau * lambda) * (6 + (6 * tau * lambda) + (3 * tau^2 * lambda^2) 
                        + (tau^3 * lambda^3)) - (p * 6)
}

wane_3months <- uniroot(Fk, c(0,1), tau = 92, p = 0.6)$root
wane_8months <- uniroot(Fk, c(0,1), tau = 244, p = 0.6)$root
# contact matrices --------------------------------------------------
path <- "/rivm/s/ainsliek/data/contact_matrices/converted/"
# path <- "inst/extdata/inputs/contact_matrices/converted/"
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
september_2020 <- readRDS(paste0(path,"transmission_matrix_september_2020.rds"))
# february_2021  <- readRDS(paste0(path,"transmission_matrix_february_2021.rds"))
# june_2021      <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))
# november_2021  <- readRDS(paste0(path,"transmission_matrix_november_2021.rds"))

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

# specify initial model parameters ---------------------------------
# parameters must be in a named list
params <- list(N = n_vec,  # population size
               # rates
               beta = 0.0004,
               beta1 = 0.14,
               sigma = 0.5,
               epsilon = 0.0,
               omega = wane_8months,
               gamma = i2r,
               h = i2h,
               i1 = h2ic,
               d = h2d,
               r = h2r,
               i2 = ic2hic,
               d_ic = ic2d,
               d_hic = hic2d,
               r_ic = hic2r,
               # simulation start time
               t_calendar_start = yday(as.Date("2020-01-01")), 
               # contact matrices for different levels of NPIs
               c_start = april_2017,
               c_lockdown = april_2020,
               c_relaxed = june_2020,
               c_very_relaxed = september_2020,
               c_normal = april_2017,
               keep_cm_fixed = FALSE,
               # IC admission thresholds
               thresh_n = 1/100000 * sum(n_vec),
               thresh_l = 3/100000 * sum(n_vec),
               thresh_m = 10/100000 * sum(n_vec),
               thresh_u = 40/100000 * sum(n_vec),
               # vaccination parameters
               eta = 1- ve_inf$mean_ve,
               eta_hosp = 1 - ve_hosp$mean_ve,
               eta_trans = 1 - ve_trans$mean_ve
              )

# Specify initial conditions --------------------------------------
empty_state <- c(rep(0, 9)) # vector of zeros

init <- c(
  t = 0,
  S = c(n_vec[1:4], n_vec[5]-1, n_vec[6:9]),
  Sv = empty_state,
  E = empty_state,
  Ev = empty_state,
  I = c(rep(0,4),1,rep(0,4)),
  Iv = empty_state,
  H = empty_state,
  Hv = empty_state,
  IC = empty_state,
  ICv = empty_state,
  H_IC = empty_state,
  H_ICv = empty_state,
  D = empty_state,
  R = empty_state,
  Rv = empty_state,
  R_1w = empty_state, 
  Rv_1w = empty_state, 
  R_2w = empty_state, 
  Rv_2w = empty_state, 
  R_3w = empty_state, 
  Rv_3w = empty_state, 
  )

# Run forward simulations --------------------------------------------
t_start <- init[1]
t_end <- t_start + 365
times <- as.integer(seq(t_start, t_end, by = 1))
betas <- readRDS("../vacamole/inst/extdata/results/model_fits/beta_draws.rds")
# sample 100 betas from last time window
betas100 <- sample(betas[[1]]$beta, 100)

# register parallel backend
registerDoParallel(cores=15)
n_sim <- 100
# Scenario 
scenarioA <- foreach(i = 1:n_sim) %dopar% {
  paramsA$beta <- betas100[i]
  paramsA$contact_mat <- april_2017[[i]]
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond, times, age_struct_seir_ode2, paramsA, method = rk45)
  as.data.frame(seir_out)
}
saveRDS(scenarioA, "/rivm/s/ainsliek/results/avacado_results.rds")

