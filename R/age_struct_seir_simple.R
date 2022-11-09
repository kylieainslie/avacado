#' Age-structured SEIR ODE model of vaccination with 2 dose primary series with 3 booster doses (3rd, 4th, and 5th doses)
#' @param times vector of times
#' @param init list of initial states
#' @param params list of parameter values
#' @return List of summary results
#' @keywords vacamole
#' @export
# Define model -----------------------------------------------------
age_struct_seir_simple <- function(times, init, params) {
  with(as.list(c(params, init)), {
    #print(t)
    # define initial state vectors from input ----------------------
    # susceptible
    S <- c(S1, S2, S3, S4, S5, S6, S7, S8, S9)
    Sv <- c(Sv1, Sv2, Sv3, Sv4, Sv5, Sv6, Sv7, Sv8, Sv9)
    # exposed
    E <- c(E1, E2, E3, E4, E5, E6, E7, E8, E9)
    Ev <- c(Ev1, Ev2, Ev3, Ev4, Ev5, Ev6, Ev7, Ev8, Ev9)
    # infectious
    I <- c(I1, I2, I3, I4, I5, I6, I7, I8, I9)
    Iv <- c(Iv1, Iv2, Iv3, Iv4, Iv5, Iv6, Iv7, Iv8, Iv9)
    # hospitalized
    H <- c(H1, H2, H3, H4, H5, H6, H7, H8, H9)
    Hv <- c(Hv1, Hv2, Hv3, Hv4, Hv5, Hv6, Hv7, Hv8, Hv9)
    # ICU
    IC <- c(IC1, IC2, IC3, IC4, IC5, IC6, IC7, IC8, IC9)
    ICv <- c(ICv1, ICv2, ICv3, ICv4, ICv5, ICv6, ICv7, ICv8, ICv9)
    # return to hospital ward after ICU
    H_IC <- c(H_IC1, H_IC2, H_IC3, H_IC4, H_IC5, H_IC6, H_IC7, H_IC8, H_IC9)
    H_ICv <- c(H_ICv1, H_ICv2, H_ICv3, H_ICv4, H_ICv5, H_ICv6, H_ICv7, H_ICv8, H_ICv9)
    # death
    D <- c(D1, D2, D3, D4, D5, D6, D7, D8, D9)
    # recovered
    R <- c(R1, R2, R3, R4, R5, R6, R7, R8, R9)
    Rv <- c(Rv1, Rv2, Rv3, Rv4, Rv5, Rv6, Rv7, Rv8, Rv9)
    # extra recovered compartments so that waning immunity is not exponential
    # R_1w <- c(R_1w1, R_1w2, R_1w3, R_1w4, R_1w5, R_1w6, R_1w7, R_1w8, R_1w9)
    # Rv_1w <- c(Rv_1w1, Rv_1w2, Rv_1w3, Rv_1w4, Rv_1w5, Rv_1w6, Rv_1w7, Rv_1w8, Rv_1w9)
    # 
    # R_2w <- c(R_2w1, R_2w2, R_2w3, R_2w4, R_2w5, R_2w6, R_2w7, R_2w8, R_2w9)
    # Rv_2w <- c(Rv_2w1, Rv_2w2, Rv_2w3, Rv_2w4, Rv_2w5, Rv_2w6, Rv_2w7, Rv_2w8, Rv_2w9)
    # 
    # R_3w <- c(R_3w1, R_3w2, R_3w3, R_3w4, R_3w5, R_3w6, R_3w7, R_3w8, R_3w9)
    # Rv_3w <- c(Rv_3w1, Rv_3w2, Rv_3w3, Rv_3w4, Rv_3w5, Rv_3w6, Rv_3w7, Rv_3w8, Rv_3w9)
    # ---------------------------------------------------------------
    
    # vaccination rate ----------------------------------------------
    if(is.null(t_vac_start) | is.null(t_vac_end)){
      alpha <- c(rep(0,9))
    } else if (!is.null(t_vac_start) & times >= t_vac_start &
               !is.null(t_vac_end) & times <= t_vac_end){
      alpha <- alpha
    } else {alpha <- c(rep(0,9))}
    
    # ---------------------------------------------------------------
    
    # determine contact matrix based on IC admissions ---------------
    ic_admissions <- sum(i1 * (H + Hv))

    # initialise flags
    if(times == 0 | params$keep_cm_fixed){
      flag_open <- 0
    }

    # determine contact matrix to use based on criteria
    tmp2 <- choose_contact_matrix_simple(params = params,
                                         criteria = ic_admissions,
                                         flag_open = flag_open,
                                         keep_fixed = params$keep_cm_fixed)
    contact_mat <- tmp2$contact_matrix
    flag_open <- tmp2$flag_open

    # determine force of infection ----------------------------------
    # seasonality
    calendar_day <- lubridate::yday(as.Date(times, origin = calendar_start_date))
    # emergence of new variants
    # var1 <- (sigma1 + sigma2*tanh((t-t_var1)/l))

    # determine transmission rate with seasonal/variant effects
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24))
      # * (var_emerg) * (1 + var1 + var2) + (!var_emerg) * 1
    # calculate force of infection (lambda)
    lambda <- beta_t * (contact_mat %*% (I + (eta_trans * Iv)))
    # ---------------------------------------------------------------
    
    # vaccination denominator
    denom <- S + R #+ R_1w + R_2w + R_3w
    #################################################################
    # ODEs:
    dS  <- (-lambda * S) - (alpha * (S/denom)) + (omega * R)
    dSv <- (alpha * (S/denom)) - (eta * lambda * Sv) + (omega * Rv)

    dE  <- (lambda * S) - (sigma * E) + epsilon
    dEv <- (eta * lambda * Sv) - (sigma * Ev)

    dI  <- (sigma * E)- (gamma * I) - (h * I)
    dIv <- (sigma * Ev) - (gamma * Iv) - (eta_hosp * h * Iv)

    dH  <- (h * I) - (i1 * H) - (d * H) - (r * H)
    dHv <- (eta_hosp * h * Iv) - (i1 * Hv) - (d * Hv) - (r * Hv)

    dIC  <- (i1 * H) - (i2 * IC) - (d_ic * IC)
    dICv <- (i1 * Hv) - (i2 * ICv) - (d_ic * ICv)

    dH_IC  <- (i2 * IC) - (r_ic * H_IC) - (d_hic * H_IC)
    dH_ICv <- (i2 * ICv) - (r_ic * H_ICv) - (d_hic * H_ICv)

    dD <- (d * (H + Hv)) + (d_ic * (IC + ICv)) + (d_hic * (H_IC + H_ICv))

    dR  <- (gamma * I) + (r * H) + (r_ic * H_IC) - alpha * (R/denom) - (omega * R) 
    dRv <- (gamma * Iv) + (r * Hv) + (r_ic * H_ICv) + alpha * (R/denom) - (omega * Rv)

    # dR_1w  <- (omega*4) * R - alpha * (R_1w/denom) - (omega*4) * R_1w
    # dRv_1w <- (omega*4) * Rv + alpha * (R_1w/denom) - (omega*4) * Rv_1w
    # 
    # dR_2w  <- (omega*4) * R_1w - alpha * (R_2w/denom) - (omega*4) * R_2w
    # dRv_2w <- (omega*4) * Rv_1w + alpha * (R_2w/denom) - (omega*4) * Rv_2w
    # 
    # dR_3w  <- (omega*4) * R_2w - alpha * (R_3w/denom) - (omega*4) * R_3w
    # dRv_3w <- (omega*4) * Rv_2w + alpha * (R_3w/denom) - (omega*4) * Rv_3w
    #################################################################
    dt <- 1
    
    # assign variables to global environment, so they can be used for
    # the next iteration
    assign("flag_open", flag_open, envir = globalenv())
    
    # output --------------------------------------------------------
    list(c(dt, dS, dSv, dE, dEv, dI, dIv, dH, dHv, dIC, dIC,
           dH_IC, dH_ICv, dD, dR, dRv#, 
           #dR_1w, dRv_1w,dR_2w, dRv_2w, dR_3w, dRv_3w
    ))
  })
}
