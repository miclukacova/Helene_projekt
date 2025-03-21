int_effect <- function(N1, N2, beta_A0_D, beta_L0_L, beta_A0_L, beta_L_D, beta_L0_D,
                       nu, eta, tau) {
  
  # Generate large data
  data0 <- sim_data_setting2(N1, 
                             eta = eta, 
                             nu = nu,
                             beta_A0_D = beta_A0_D, 
                             beta_L0_L = beta_L0_L, 
                             beta_A0_L = beta_A0_L, 
                             beta_L_D = beta_L_D, 
                             beta_L0_D = beta_L0_D,
                             cens = 0)
  
  # Group data based on treatment
  no_treat_group <- data0[A0 == 0]
  treat_group <- data0[A0 == 1]
  
  # Fit Weibull
  survfit <- survreg(Surv(Time, Delta == 3) ~ 1, 
                     data = no_treat_group[L == 0], 
                     cluster = ID,
                     dist='weibull')
  
  # Estimates in no treatment group
  nu_est <- 1/survfit$scale
  eta_est <- 1/(exp(survfit$coefficients[1]))^nu_est
  
  # Generate large data set under the intervened intensity 
  data_new <- sim_data_setting2(N = N2, 
                                cens = 0,
                                eta = c(eta[1:3],eta_est),
                                nu = c(nu[1:3],nu_est),
                                beta_A0_D = beta_A0_D, 
                                beta_L_D = beta_L_D, 
                                beta_L0_D = beta_L0_D,
                                beta_L0_L = 0,
                                beta_A0_L = 0)

  #Generate large data set without intervened intensity
  data_new_no_int <- sim_data_setting2(N = N2, 
                                       eta = eta, 
                                       nu = nu,
                                       beta_A0_D = beta_A0_D, 
                                       beta_L0_L = beta_L0_L, 
                                       beta_A0_L = beta_A0_L, 
                                       beta_L_D = beta_L_D, 
                                       beta_L0_D = beta_L0_D,
                                       cens = 0)
  #Proportion of subjects dying before some time $\tau$ in treatment group
  prop_int <- mean(data_new[Delta == 1 & A0 == 1, Time] < tau) # with intervention
  prop_no_int <- mean(data_new_no_int[Delta == 1 & A0 == 1, Time] < tau) # without intervention
  
  return(c(prop_int, prop_no_int))
}



compare_int_effect <- function(N1, N2, beta_A0_D = 0, beta_L0_L, beta_A0_L, 
                               beta_L_D = 0, beta_L0_D,
                               nu, eta, tau, setting = 1) {
  
  # We find out the length of the varying parameter
  B <-  max(length(beta_L0_L), length(beta_L_D), length(beta_A0_L), length(beta_L0_D), length(beta_A0_D))
  
  # All other parameters are repeated 
  if(length(beta_L0_L) < B) beta_L0_L <- rep(beta_L0_L, B)
  if(length(beta_L_D) < B) beta_L_D <- rep(beta_L_D, B)
  if(length(beta_A0_L) < B) beta_A0_L <- rep(beta_A0_L, B)
  if(length(beta_L0_D) < B) beta_L0_D <- rep(beta_L0_D, B)
  if(length(beta_A0_D) < B) beta_A0_D <- rep(beta_A0_D, B)
  
  # A matrix for the results
  res <- matrix(nrow = B, ncol = 2*3)
  
  if(setting == 1){
    for(i in 1:B){
      res[i,1:2] <- int_effect(N1 = N1, 
                               N2 = N2,
                               eta = eta, 
                               nu = nu,
                               beta_A0_D = 0, 
                               beta_L0_L = beta_L0_L[i], 
                               beta_A0_L = beta_A0_L[i], 
                               beta_L_D = beta_L_D[i], 
                               beta_L0_D = beta_L0_D[i],
                               tau = tau)
      
      res[i,3:4] <- int_effect(N1 = N1, 
                               N2 = N2,
                               eta = eta, 
                               nu = nu,
                               beta_A0_D = -0.1, 
                               beta_L0_L = beta_L0_L[i], 
                               beta_A0_L = beta_A0_L[i], 
                               beta_L_D = beta_L_D[i], 
                               beta_L0_D = beta_L0_D[i],
                               tau = tau)
      
      res[i,5:6] <- int_effect(N1 = N1, 
                               N2 = N2,
                               eta = eta, 
                               nu = nu,
                               beta_A0_D = -0.2, 
                               beta_L0_L = beta_L0_L[i], 
                               beta_A0_L = beta_A0_L[i], 
                               beta_L_D = beta_L_D[i], 
                               beta_L0_D = beta_L0_D[i],
                               tau = tau)
    }
  }
  else{
    for(i in 1:B){
      res[i,1:2] <- int_effect(N1 = N1, 
                                 N2 = N2,
                                 eta = eta, 
                                 nu = nu,
                                 beta_A0_D = beta_A0_D[i], 
                                 beta_L0_L = beta_L0_L[i], 
                                 beta_A0_L = beta_A0_L[i], 
                                 beta_L_D = 0, 
                                 beta_L0_D = beta_L0_D[i],
                                 tau = tau)
      res[i,3:4] <- int_effect(N1 = N1, 
                                 N2 = N2,
                                 eta = eta, 
                                 nu = nu,
                                 beta_A0_D = beta_A0_D[i], 
                                 beta_L0_L = beta_L0_L[i], 
                                 beta_A0_L = beta_A0_L[i], 
                                 beta_L_D = 0.5, 
                                 beta_L0_D = beta_L0_D[i],
                                 tau = tau)
      res[i,5:6] <- int_effect(N1 = N1, 
                                 N2 = N2,
                                 eta = eta, 
                                 nu = nu,
                                 beta_A0_D = beta_A0_D[i], 
                                 beta_L0_L = beta_L0_L[i], 
                                 beta_A0_L = beta_A0_L[i], 
                                 beta_L_D = 1, 
                                 beta_L0_D = beta_L0_D[i],
                                 tau = tau)
      }
    }
  return(res)
}