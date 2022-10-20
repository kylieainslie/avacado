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
    R_1w <- c(R_1w1, R_1w2, R_1w3, R_1w4, R_1w5, R_1w6, R_1w7, R_1w8, R_1w9)
    Rv_1w <- c(Rv_1w1, Rv_1w2, Rv_1w3, Rv_1w4, Rv_1w5, Rv_1w6, Rv_1w7, Rv_1w8, Rv_1w9)

    R_2w <- c(R_2w1, R_2w2, R_2w3, R_2w4, R_2w5, R_2w6, R_2w7, R_2w8, R_2w9)
    Rv_2w <- c(Rv_2w1, Rv_2w2, Rv_2w3, Rv_2w4, Rv_2w5, Rv_2w6, Rv_2w7, Rv_2w8, Rv_2w9)

    R_3w <- c(R_3w1, R_3w2, R_3w3, R_3w4, R_3w5, R_3w6, R_3w7, R_3w8, R_3w9)
    Rv_3w <- c(Rv_3w1, Rv_3w2, Rv_3w3, Rv_3w4, Rv_3w5, Rv_3w6, Rv_3w7, Rv_3w8, Rv_3w9)

    # define vaccination parameters ---------------------------------
    # don't index parameters when there's no vaccination, it's faster!
    index <- floor(times) + 1              # use floor of time point + 1 to index df
    # daily vac rate
    alpha1 <- params$alpha1[index, -1] # remove date column
    # delay to protection
    delay1 <- params$delay1[index, -1]
    # protection against infection (1 - VE_inf)
    eta1   <- params$eta1[index, -1]
    # protection against hospitalisation 1 - (1 - VE_hosp) / (1 - VE_inf)
    eta_hosp1   <- params$eta_hosp1[index, -1]
    # protection against transmission (1 - VE_trans)
    eta_trans1   <- as.numeric(params$eta_trans1[index, -1])
    # ---------------------------------------------------------------

    # determine force of infection ----------------------------------
    # seasonality
    calendar_day <- lubridate::yday(as.Date(times, origin = calendar_start_date))
    # emergence of new variants
    var1 <- (sigma1 + sigma2*tanh((t-t_var1)/l))
    var2 <- (sigma1 + sigma2*tanh((t-t_var2)/l))
    # determine transmission rate with seasonal/variant effects
    beta_t <- beta * (1 + beta1 * cos(2 * pi * calendar_day / 365.24)) *
      (var_emerg) * (1 + var1 + var2) + (!var_emerg) * 1
    # calculate force of infection (lambda)
    lambda <- beta_t * (contact_mat %*% (I + (eta_trans1 * Iv_1d) + (eta_trans2 * Iv_2d) + (eta_trans3 * Iv_3d)
                                         + (eta_trans4 * Iv_4d) + (eta_trans5 * Iv_5d)
    ))
    # ---------------------------------------------------------------

    #################################################################
    # ODEs:
    dS <- -lambda * S - alpha1 * S + (omega*4) * R_3w
    dShold_1d <- alpha1 * S - (1/delay1) * Shold_1d - lambda * Shold_1d
    dSv_1d <- (1/delay1) * Shold_1d - eta1 * lambda * Sv_1d - alpha2 * Sv_1d + (omega*4) * Rv_1d_3w
    dShold_2d <- alpha2 * Sv_1d - (1/delay2) * Shold_2d - eta1 * lambda * Shold_2d
    dSv_2d <- (1/delay2) * Shold_2d - eta2 * lambda * Sv_2d - alpha3 * Sv_2d + (omega*4) * Rv_2d_3w
    dShold_3d <- alpha3 * Sv_2d - (1/delay3) * Shold_3d - eta2 * lambda * Shold_3d
    dSv_3d <- (1/delay3) * Shold_3d - eta3 * lambda * Sv_3d - alpha4 * Sv_3d + (omega*4) * Rv_3d_3w
    dShold_4d <- alpha4 * Sv_3d - (1/delay4) * Shold_4d - eta3 * lambda * Shold_4d
    dSv_4d <- (1/delay4) * Shold_4d - eta4 * lambda * Sv_4d - alpha5 * Sv_4d + (omega*4) * Rv_4d_3w
    dShold_5d <- alpha5 * Sv_4d - (1/delay5) * Shold_5d - eta4 * lambda * Shold_5d
    dSv_5d <- (1/delay5) * Shold_5d - eta5 * lambda * Sv_5d + (omega*4) * Rv_5d_3w

    dE     <- lambda * S + lambda * Shold_1d - sigma * E + epsilon
    dEv_1d <- eta1 * lambda * Sv_1d + eta1 * lambda * Shold_2d - sigma * Ev_1d
    dEv_2d <- eta2 * lambda * Sv_2d + eta2 * lambda * Shold_3d - sigma * Ev_2d
    dEv_3d <- eta3 * lambda * Sv_3d + eta3 * lambda * Shold_4d - sigma * Ev_3d
    dEv_4d <- eta4 * lambda * Sv_4d + eta4 * lambda * Shold_5d - sigma * Ev_4d
    dEv_5d <- eta5 * lambda * Sv_5d - sigma * Ev_5d

    dI     <- sigma * E - (gamma + h) * I
    dIv_1d <- sigma * Ev_1d - (gamma + eta_hosp1 * h) * Iv_1d
    dIv_2d <- sigma * Ev_2d - (gamma + eta_hosp2 * h) * Iv_2d
    dIv_3d <- sigma * Ev_3d - (gamma + eta_hosp3 * h) * Iv_3d
    dIv_4d <- sigma * Ev_4d - (gamma + eta_hosp4 * h) * Iv_4d
    dIv_5d <- sigma * Ev_5d - (gamma + eta_hosp5 * h) * Iv_5d

    dH     <- h * I - (i1 + d + r) * H
    dHv_1d <- eta_hosp1 * h * Iv_1d - (i1 + d + r) * Hv_1d
    dHv_2d <- eta_hosp2 * h * Iv_2d - (i1 + d + r) * Hv_2d
    dHv_3d <- eta_hosp3 * h * Iv_3d - (i1 + d + r) * Hv_3d
    dHv_4d <- eta_hosp4 * h * Iv_4d - (i1 + d + r) * Hv_4d
    dHv_5d <- eta_hosp5 * h * Iv_5d - (i1 + d + r) * Hv_5d

    dIC     <- i1 * H - (i2 + d_ic) * IC
    dICv_1d <- i1 * Hv_1d - (i2 + d_ic) * ICv_1d
    dICv_2d <- i1 * Hv_2d - (i2 + d_ic) * ICv_2d
    dICv_3d <- i1 * Hv_3d - (i2 + d_ic) * ICv_3d
    dICv_4d <- i1 * Hv_4d - (i2 + d_ic) * ICv_4d
    dICv_5d <- i1 * Hv_5d - (i2 + d_ic) * ICv_5d

    dH_IC     <- i2 * IC - (r_ic + d_hic) * H_IC
    dH_ICv_1d <- i2 * ICv_1d - (r_ic + d_hic) * H_ICv_1d
    dH_ICv_2d <- i2 * ICv_2d - (r_ic + d_hic) * H_ICv_2d
    dH_ICv_3d <- i2 * ICv_3d - (r_ic + d_hic) * H_ICv_3d
    dH_ICv_4d <- i2 * ICv_4d - (r_ic + d_hic) * H_ICv_4d
    dH_ICv_5d <- i2 * ICv_5d - (r_ic + d_hic) * H_ICv_5d

    dD <- d * (H + Hv_1d + Hv_2d + Hv_3d + Hv_4d + Hv_5d) +               #
      d_ic * (IC + ICv_1d + ICv_2d + ICv_3d + ICv_4d + ICv_5d) +           #
      d_hic * (H_IC + H_ICv_1d + H_ICv_2d + H_ICv_3d + H_ICv_4d + H_ICv_5d) #

    dR     <- (gamma * I) + (r * H) + (r_ic * H_IC) - ((omega*4) * R)
    dRv_1d <- (gamma * Iv_1d) + (r * Hv_1d) + (r_ic * H_ICv_1d) - ((omega*4) * Rv_1d)
    dRv_2d <- (gamma * Iv_2d) + (r * Hv_2d) + (r_ic * H_ICv_2d) - ((omega*4) * Rv_2d)
    dRv_3d <- (gamma * Iv_3d) + (r * Hv_3d) + (r_ic * H_ICv_3d) - ((omega*4) * Rv_3d)
    dRv_4d <- (gamma * Iv_4d) + (r * Hv_4d) + (r_ic * H_ICv_4d) - ((omega*4) * Rv_4d)
    dRv_5d <- (gamma * Iv_5d) + (r * Hv_5d) + (r_ic * H_ICv_5d) - ((omega*4) * Rv_5d)

    dR_1w     <- (omega*4) * R - (omega*4) * R_1w
    dRv_1d_1w <- (omega*4) * Rv_1d - (omega*4) * Rv_1d_1w
    dRv_2d_1w <- (omega*4) * Rv_2d - (omega*4) * Rv_2d_1w
    dRv_3d_1w <- (omega*4) * Rv_3d - (omega*4) * Rv_3d_1w
    dRv_4d_1w <- (omega*4) * Rv_4d - (omega*4) * Rv_4d_1w
    dRv_5d_1w <- (omega*4) * Rv_5d - (omega*4) * Rv_5d_1w

    dR_2w     <- (omega*4) * R_1w - (omega*4) * R_2w
    dRv_1d_2w <- (omega*4) * Rv_1d_1w - (omega*4) * Rv_1d_2w
    dRv_2d_2w <- (omega*4) * Rv_2d_1w - (omega*4) * Rv_2d_2w
    dRv_3d_2w <- (omega*4) * Rv_3d_1w - (omega*4) * Rv_3d_2w
    dRv_4d_2w <- (omega*4) * Rv_4d_1w - (omega*4) * Rv_4d_2w
    dRv_5d_2w <- (omega*4) * Rv_5d_1w - (omega*4) * Rv_5d_2w

    dR_3w     <- (omega*4) * R_2w - (omega*4) * R_3w
    dRv_1d_3w <- (omega*4) * Rv_1d_2w - (omega*4) * Rv_1d_3w
    dRv_2d_3w <- (omega*4) * Rv_2d_2w - (omega*4) * Rv_2d_3w
    dRv_3d_3w <- (omega*4) * Rv_3d_2w - (omega*4) * Rv_3d_3w
    dRv_4d_3w <- (omega*4) * Rv_4d_2w - (omega*4) * Rv_4d_3w
    dRv_5d_3w <- (omega*4) * Rv_5d_2w - (omega*4) * Rv_5d_3w
    #################################################################
    dt <- 1
    # output --------------------------------------------------------
    list(c(dt, dS, dSv_1d, dSv_2d,dSv_3d, dSv_4d, dSv_5d,
           dShold_1d, dShold_2d,  dShold_3d, dShold_4d, dShold_5d,
           dE, dEv_1d, dEv_2d, dEv_3d, dEv_4d, dEv_5d,
           dI, dIv_1d, dIv_2d, dIv_3d, dIv_4d, dIv_5d,
           dH, dHv_1d, dHv_2d, dHv_3d, dHv_4d, dHv_5d,
           dIC, dICv_1d, dICv_2d, dICv_3d, dICv_4d, dICv_5d,
           dH_IC, dH_ICv_1d, dH_ICv_2d, dH_ICv_3d, dH_ICv_4d, dH_ICv_5d,
           dD,
           dR, dRv_1d, dRv_2d, dRv_3d, dRv_4d, dRv_5d,
           dR_1w, dRv_1d_1w, dRv_2d_1w, dRv_3d_1w, dRv_4d_1w, dRv_5d_1w,
           dR_2w, dRv_1d_2w, dRv_2d_2w, dRv_3d_2w, dRv_4d_2w, dRv_5d_2w,
           dR_3w, dRv_1d_3w, dRv_2d_3w, dRv_3d_3w, dRv_4d_3w, dRv_5d_3w
    ))
  })
}
