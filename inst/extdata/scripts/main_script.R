# Script for running scenarios for the Eupropean Scenario Hub
# URL: https://github.com/covid19-forecast-hub-europe/covid19-scenario-hub-europe#readme

# preamble ---------------------------------------------------------
# This script will load necessary packages, data sources, fit the model to data,
# and then run scenarios.
# All scenarios will be using Dutch data
# TODO: generalize to other European countries
# TODO: break up code into work chunks
# ------------------------------------------------------------------

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

devtools::load_all() # load vacamole
library(vacamole)
# -------------------------------------------------------------------

# Load data ---------------------------------------------------------
# if off the server, read in from inst/extdata/data
data_date <- "2022-05-22"
osiris1 <- readRDS(paste0("inst/extdata/data/case_data_upto_", data_date, ".rds"))

# read in transition rates -----------------------------------------
transition_rates <- readRDS("inst/extdata/inputs/transition_rates.rds")

# define population size (by age group)
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 
              0.13083463,0.14514332, 0.12092904, 0.08807406, 
              0.04622194)
n <- 17407585 # Dutch population size
n_vec <- n * age_dist

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
#path <- "/rivm/s/ainsliek/data/contact_matrices/converted/"
path <- "inst/extdata/inputs/contact_matrices/converted/"
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
september_2020 <- readRDS(paste0(path,"transmission_matrix_september_2020.rds"))
february_2021  <- readRDS(paste0(path,"transmission_matrix_february_2021.rds"))
june_2021      <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))
november_2021  <- readRDS(paste0(path,"transmission_matrix_november_2021.rds"))

# put contact matrices into a list for input into fit_to_data_func()
cm_list <- list(
  april_2017 = april_2017,
  april_2020 = april_2020,
  june_2020 = june_2020,
  september_2020 = september_2020,
  february_2021 = february_2021,
  june_2021 = june_2021,
  november_2021 = november_2021
)
# ve estimates ------------------------------------------------------
# list containing the following named lists:
# delays, ve_inf, ve_hosp, ve_trans
# each named list has the following named elements:
# pfizer, moderna, astrazeneca, jansen
ve_params <- readRDS("inst/extdata/inputs/ve_params.rds")

# specify initial model parameters ---------------------------------
# parameters must be in a named list
params <- list(beta = 0.0004,
               beta1 = 0.14,
               gamma = 0.5,
               sigma = 0.5,
               epsilon = 0.0,
               omega = wane_8months,
               N = n_vec,
               h = transition_rates$h,
               i1 = transition_rates$i1,
               i2 = transition_rates$i2,
               d = transition_rates$d, 
               d_ic = transition_rates$d_ic,
               d_hic = transition_rates$d_hic,
               r = transition_rates$r,
               r_ic = transition_rates$r_ic,
               p_report = 1/3,
               c_start = april_2017,
               keep_cm_fixed = FALSE,
               #vac_inputs = vac_rates_wt,
               use_cases = FALSE,  
               no_vac = FALSE,
               t_calendar_start = yday(as.Date("2020-01-01")), 
               beta_change = 0.0001,
               t_beta_change = 165
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
times <- seq(0,250, by = 1)

seir_out <- lsoda(init, times, age_struct_seir_ode2, params)
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output2(seir_out)

