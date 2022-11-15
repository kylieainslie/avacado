#' Summarise output from SEIR model for plotting
#' @param seir_output list of dataframes from post-processing
#' @param params list of parameter values
#' @param start_date calendar date of start of simulation
#' @param times vector of time points
#' @param vac_inputs vaccination inputs (from the output of convert_vac_schedule())
#' @return List of summary results
#' @keywords vacamole
#' @export

summarise_results_simple <- function(seir_output, params, t_vec) {
  
  # determine contact matrix for each time step (to calculate correct FOI below)
  contact_mat_list <- list()
  # determine contact matrix based on IC admissions
  ic_admin <- rowSums(params$i1 * (seir_output$H + seir_output$Hv))
  
  for(t in 1:length(t_vec)){
    # initialise flags
    if(t == 1 | params$keep_cm_fixed){
      flag_open <- 0
    }
    
    # determine contact matrix to use based on criteria
    tmp2 <- choose_contact_matrix_simple(params = params, 
                                         criteria = ic_admin[t], 
                                         flag_open = flag_open, 
                                         keep_fixed = params$keep_cm_fixed)
    contact_mat_list[[t]] <- tmp2$contact_matrix
    flag_open <- tmp2$flag_open
    
  }
  
  # get force of infection (lambda) --------------------------------------------
  calendar_day <- lubridate::yday(as.Date(t_vec, 
                                          origin = params$calendar_start_date))
  beta_t <- params$beta * (1 + params$beta1 * cos(2 * pi * calendar_day / 365.24)) 
  lambda <- get_foi_simple(x  = seir_output, 
                           y1 = params$eta_trans,
                           beta = beta_t, 
                           contact_mat = contact_mat_list,
                           times = t_vec)
  
  # calculate infections -------------------------------------------------------
  new_infections <- lambda * (seir_output$S + (params$eta * (seir_output$Sv))) %>%
    rename_with(., ~ paste0("age_group",1:9))
  
  new_infections <- new_infections %>%
    mutate(target_variable = "inc infection",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  # calculate cases ------------------------------------------------------------
  # new_cases <- sweep(params$sigma * (seir_output$E + seir_output$Ev_1d + 
  #   seir_output$Ev_2d + seir_output$Ev_3d + seir_output$Ev_4d + seir_output$Ev_5d),
  #   2, params$p_report, "*") %>%
  #   rename_with(., ~ paste0("age_group",1:9)) 
  # 
  # new_cases <- new_cases %>%
  #   mutate(target_variable = "inc case",
  #          time = t_vec,
  #          date = as.Date(time, origin = params$calendar_start_date),
  #          epiweek = lubridate::epiweek(date),
  #          year = lubridate::epiyear(date),
  #          wk = floor(difftime(date, date[1], units = "weeks")) + 1,
  #          horizon = paste(wk, "wk")) %>%
  #   select(-wk)
  # 
  # calculate hospital admissions ----------------------------------------------
  hosp_admissions <- (sweep(seir_output$I, 2, params$h, "*") +
    sweep(seir_output$Iv, 2, params$h_v, "*")) %>%
    rename_with(., ~ paste0("age_group",1:9)) 
  
  hosp_admissions <- hosp_admissions %>%
    mutate(target_variable = "inc hosp",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # calculate hospital occupancy -----------------------------------------------
  hosp_occ <-  (seir_output$H + seir_output$H_IC + seir_output$Hv + seir_output$H_ICv) %>%
    rename_with(., ~ paste0("age_group",1:9)) 
  
  hosp_occ <- hosp_occ %>%
    mutate(target_variable = "occ hosp",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # calculate IC admissions ----------------------------------------------------
  ic_admissions <- sweep(seir_output$H + seir_output$Hv, 2, params$i1, "*") %>%
    rename_with(., ~ paste0("age_group",1:9)) 
  
  ic_admissions <- ic_admissions %>%
    mutate(target_variable = "inc icu",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # calculate hospital occupancy -----------------------------------------------
  ic_occ <-  (seir_output$IC + seir_output$ICv) %>%
    rename_with(., ~ paste0("age_group",1:9)) 
  
  ic_occ <- ic_occ %>%
    mutate(target_variable = "occ icu",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # calculate deaths -----------------------------------------------------------
  new_deaths <- 
    # from hospital compartment
    sweep(seir_output$H + seir_output$Hv, 2, params$d, "*") +
    # from IC compartment
    sweep(seir_output$IC + seir_output$ICv, 2, params$d_ic, "*") + 
    # from hospital after IC
    sweep(seir_output$H_IC + seir_output$H_ICv, 2, params$d_hic, "*")
   
  
  new_deaths <- new_deaths %>%
    rename_with(., ~ paste0("age_group",1:9)) %>%
    mutate(target_variable = "inc death",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # calculate number immune (people in R compartment)
  recovered <-  (seir_output$R + seir_output$Rv) %>%
    rename_with(., ~ paste0("age_group",1:9)) 
  
  recovered <- recovered %>%
    mutate(target_variable = "occ rec",
           time = t_vec,
           date = as.Date(time, origin = params$calendar_start_date),
           epiweek = lubridate::epiweek(date),
           year = lubridate::epiyear(date),
           wk = floor(difftime(date, date[1], units = "weeks")) + 1,
           horizon = paste(wk, "wk")) %>%
    select(-wk)
  
  # Create object into format for scenario hub ---------------------------------
  # bind all outcome data frames together and wrangle
  rtn <- bind_rows(new_infections, 
                   # new_cases, 
                   hosp_admissions, 
                   hosp_occ,
                   ic_admissions, 
                   ic_occ,
                   new_deaths,
                   recovered) %>%
    pivot_longer(cols = age_group1:age_group9,
                 names_to = "age_group",
                 names_prefix = "age_group",
                 values_to = "value") %>%
    select(time, date, epiweek, year, horizon, target_variable, age_group, value)

  return(rtn)
}
