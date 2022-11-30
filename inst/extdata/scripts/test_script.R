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
               beta = 0.0004,
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

# Specify initial conditions -----------------------------------------
# get initial conditions from last time point (275) of Wave 1
# read in results from Wave 1
scenarioB_wave1 <- readRDS("/rivm/s/ainsliek/results/covid_scenarios/wave1_scenarioB.rds")

# get tail of dataset for each simulation
initB <- lapply(scenarioB_wave1, tail, 1)

# Run forward simulations --------------------------------------------
# Scenario B: voluntary
t_start <- 0 #as.numeric(initA[[1]][1])
t_end <- t_start + 365
times <- as.integer(seq(t_start, t_end, by = 1))
n_sim <- 10
# Scenario B: voluntary ----
registerDoParallel(cores=15)
# cl <- makeCluster(15)
# registerDoParallel(cl)
# clusterExport(cl,list("flag_open","rkMethod", "ode"),envir=globalenv())

scenarioB <- list()
# scenarioB <- foreach(i = 1:n_sim) %dopar% {
for (i in 1:n_sim){
  flag_open <- 0
  init <- unlist(initB[[i]][-1])
  init[1] <- t_start
  params$keep_cm_fixed <- FALSE
  #params$beta <- betas100[i] * 1.5
  params$c_start <- june_2022[[i]]
  params$c_lockdown <- april_2020[[i]]
  params$c_open <- params$c_lockdown
  
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, times, age_struct_seir_simple, params, method = rk45)
  scenarioB[[i]] <- as.data.frame(seir_out)
}
# doParallel::stopImplicitCluster()

# Post-process raw output ---------------------------------------------
outB <- list()
# loop over samples and summarize results for each scenario
for(s in 1:n_sim){
  # specify shared parameter values (transmission rate and starting contact matrix)
  #params$beta <- betas100[s] * 1.5
  params$c_start <- june_2022[[s]] # Pico 8 when available
  
  # Scenario B - voluntary
  params$keep_cm_fixed <- FALSE
  params$c_lockdown <- june_2020[[s]]
  params$c_open <- params$c_lockdown
  
  seir_outputB <- postprocess_age_struct_model_output_simple(scenarioB[[s]])
  seir_outcomesB <- summarise_results_simple(seir_outputB, params = params, t_vec = times) %>%
    mutate(sample = s)
  outB[[s]] <- seir_outcomesB
}

# data frames
dfB <- bind_rows(outB) %>% mutate(scenario_id = "B", wave = "Wave 2") 

# Plot -----------------------------------------------------------------
df_sum_ag <- dfB %>%
  group_by(wave, scenario_id, target_variable, date, sample) %>%
  summarise(sum = sum(value)) %>%
  ungroup()

# plot
p_lines <- ggplot(data = df_sum_ag %>% filter(sample %in% c(1,3,7,8,9),
                                              target_variable == "inc infection"), 
                  aes(x = date, y = sum, group = sample,
                                     color = scenario_id)) +
  geom_line() +
  facet_grid(sample~wave)

p_lines

# Average with confidence bounds -----------------------------------------------
df_summary <- df_sum_ag %>%
  group_by(wave, scenario_id, target_variable, date) %>%
  summarise(median  = median(sum),
            q025 = quantile(sum, probs = 0.025),
            q25  = quantile(sum, probs = 0.25),
            q75  = quantile(sum, probs = 0.75),
            q975 = quantile(sum, probs = 0.975)
  ) %>%
  select(date, wave, scenario_id, target_variable, median:q975) %>%
  mutate(target_variable = factor(target_variable))

# plot 
p_ribbon <- ggplot(data = df_summary %>%
                     filter(target_variable %in% c(#"inc icu",
                       #"inc infection"
                       "occ hosp",
                       "occ icu", 
                       "occ rec"
                     ),
                     wave == "Wave 2"#,
                     #scenario_id != "A"
                     ),
                   aes(x = date, y = median, color = scenario_id, fill = scenario_id)) +
  geom_line() +
  #geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.1, color = NA) +
  xlab("Date") + 
  ylab("Median value") +
  scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14,face="bold")) +
  guides(fill=guide_legend("Scenario ID"), colour = guide_legend("Scenario ID")) +
  facet_grid(target_variable~., scales = "free") #+
p_ribbon


