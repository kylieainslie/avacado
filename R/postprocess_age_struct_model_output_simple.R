#TODO refer to manuscript or somewhere else for explanation of each compartment symbol
#' Postprocess output from age-structured SEIR ODE model of vaccination
#' with 2 doses and delay to protection
#' @param dat output from seir model as a data frame
#' @return List of summary results
#' @keywords vacamole
#' @importFrom stringr str_detect
#' @export
postprocess_age_struct_model_output_simple <- function(dat) {
  # S ----------
  S <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "S[:digit:]")])
  Sv <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Sv[:digit:]")])
  
  # E ----------
  E <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "E[:digit:]")])
  Ev <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Ev[:digit:]")])
  
  # I ----------
  I <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "I[:digit:]")])
  Iv <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Iv[:digit:]")])
  
  # H ----------
  H <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H[:digit:]")])
  Hv <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Hv[:digit:]")])
  
  # H_IC -------
  H_IC <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_IC[:digit:]")])
  H_ICv <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "H_ICv[:digit:]")])
  
  # IC ---------
  IC <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "IC[:digit:]")]) %>%
    select(-names(dat)[str_detect(names(dat), pattern = "H_IC[:digit:]")])
  ICv <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "ICv[:digit:]")]) %>%
    select(-names(dat)[str_detect(names(dat), pattern = "H_ICv[:digit:]")])
  
  # D ----------
  D <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "D[:digit:]")])
  
  # R ----------
  R <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R[:digit:]")])
  Rv <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv[:digit:]")])
  
  # R_1w ------
  R_1w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R_1w[:digit:]")])
  Rv_1w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_1w[:digit:]")])
  
  # R_2w ------
  R_2w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R_2w[:digit:]")])
  Rv_2w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_2w[:digit:]")])
  
  # R_3w ------
  R_3w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "R_3w[:digit:]")])
  Rv_3w <- dat %>%
    select(names(dat)[str_detect(names(dat), pattern = "Rv_3w[:digit:]")])
  
  # output ----
  rtn <- list(
    S = S,
    Sv = Sv,
    E = E,
    Ev = Ev,
    I = I,
    Iv = Iv,
    H = H,
    Hv = Hv,
    H_IC = H_IC,
    H_ICv = H_ICv,
    IC = IC,
    ICv = ICv,
    D = D,
    R = R,
    Rv = Rv,
    R_1w = R_1w,
    Rv_1w = Rv_1w,
    R_2w = R_2w,
    Rv_2w = Rv_2w,
    R_3w = R_3w,
    Rv_3w = Rv_3w
  )

  return(rtn)
}
