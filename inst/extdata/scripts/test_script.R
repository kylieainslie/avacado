times <- seq(0, l00, by = 1)

params <- list(N = n_vec,  # population size
               # rates
               beta = 0.0004,
               beta1 = 0.14,
               sigma = 0.5,
               epsilon = 0.00,
               omega = 0.0038,
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
               c_start = april_2017$mean,
               c_lockdown = april_2020$mean,
               c_open = april_2020$mean,
               keep_cm_fixed = FALSE,
               # IC admission thresholds
               thresh_o = 1,
               thresh_l = 40,
               # vaccination parameters
               alpha = c(rep(0,9)),
               t_vac_start = NULL,
               t_vac_end = NULL,
               eta = 1, # - ve_inf$mean_ve,
               #eta_hosp = 1, # - ve_hosp$mean_ve,
               eta_trans = 1 # - ve_trans$mean_ve
)

rk45 <- rkMethod("rk45f")
seir_out <- ode(init, times, age_struct_seir_simple, params, method = rk45)
test <- as.data.frame(seir_out)

# output error message if sum of all compartments is not equal to the total population size
  if(!isTRUE(all.equal(sum(test[dim(test)[1],-c(1,2,138:146)]),sum(params$N)))){
    stop("Error: sum of compartments is not equal to population size")
  }
  
test_s <- test %>%
  select(starts_with("S"))
# head(test_s)

plot(1:nrow(test_s), rowSums(test_s), type = "l")

test_i <- test %>%
  select(starts_with("I"))
# head(test_i)

plot(1:nrow(test_i), rowSums(test_i), type = "l")
