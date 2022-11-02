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
  
  # get force of infection (lambda) --------------------------------------------
  calendar_day <- lubridate::yday(as.Date(t_vec, origin = params$calendar_start_date))
  beta_t <- params$beta * (1 + params$beta1 * cos(2 * pi * calendar_day / 365.24)) 
  lambda <- get_foi_simple(x  = seir_output, 
                           y1 = params$eta_trans,
                           beta = beta_t, 
                           contact_mat = params$c_start,
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
  hosp_admissions <- sweep(seir_output$I + params$eta_hosp * seir_output$Iv, 
                           2, params$h, "*") %>%
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
  
  # Create object into format for scenario hub ---------------------------------
  # bind all outcome data frames together and wrangle
  rtn <- bind_rows(new_infections, 
                   # new_cases, 
                   hosp_admissions, 
                   ic_admissions, 
                   new_deaths) %>%
    pivot_longer(cols = age_group1:age_group9,
                 names_to = "age_group",
                 names_prefix = "age_group",
                 values_to = "value") %>%
    select(time, date, epiweek, year, horizon, target_variable, age_group, value)

  return(rtn)
}
