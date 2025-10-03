# Time varying effects and tests

library(survival)
library(simevent)

# Tjek af at det virker når vi ikke har time varying effects

N <- 5000

eta <- rep(0.1, 4)
term_deltas <- c(0,1)
time_var_eff <- matrix(0.1, ncol = 4, nrow = 6)
#beta <- matrix(nrow = 6, ncol = 4, 0.1)
beta <- matrix(nrow = 6, ncol = 4, 10)
# t_prime bliver sat højt så vi altid er i det sædvanlige scenarie
t_prime <- 1000                                                

data <- simEventTV(N = N, t_prime = t_prime, time_var_eff = time_var_eff, eta = eta,
                            term_deltas = term_deltas, beta = beta, lower = 10^(-20), upper = 10^2,
                            max_events = 5)

data <- IntFormatData(data)

survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + N2 + N3, data = data)
survfit1 <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + N2 + N3, data = data)
survfit2 <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A0 + N2 + N3, data = data)
survfit3 <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + A0 + N2 + N3, data = data)

any(abs(survfit0$coefficients - beta[c(1,2,5,6),1]) > 0.1)
any(abs(survfit1$coefficients - beta[c(1,2,5,6),2]) > 0.1)
any(abs(survfit2$coefficients - beta[c(1,2,5,6),3]) > 0.1)
any(abs(survfit3$coefficients - beta[c(1,2,5,6),4]) > 0.1)

# Effekterne bliver sat til 0 så vi igen er i det sævanlige scenarie
t_prime <- 1
time_var_eff <- matrix(0, ncol = 4, nrow = 6)
data <- simEventTV(N = N, t_prime = t_prime, time_var_eff = time_var_eff, eta = eta,
                            term_deltas = term_deltas, beta = beta, lower = 10^(-20), upper = 10^2,
                            max_events = 5)

data <- IntFormatData(data)
survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + N2 + N3, data = data)
survfit1 <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + N2 + N3, data = data)
survfit2 <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A0 + N2 + N3, data = data)
survfit3 <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + A0 + N2 + N3, data = data)

any(abs(survfit0$coefficients - beta[c(1,2,5,6),1]) > 0.1)
any(abs(survfit1$coefficients - beta[c(1,2,5,6),2]) > 0.1)
any(abs(survfit2$coefficients - beta[c(1,2,5,6),3]) > 0.1)
any(abs(survfit3$coefficients - beta[c(1,2,5,6),4]) > 0.1)


# Tjek af at det virker når vi har time varying effects
time_var_eff <- matrix(0.3, ncol = 4, nrow = 6)
t_prime <- 1

data <- simEventTV(N = N, t_prime = t_prime, time_var_eff = time_var_eff, eta = eta,
                            term_deltas = term_deltas, beta = beta, lower = 10^(-20), upper = 10^2,
                            max_events = 5)

data <- IntFormatData(data, timeVar = TRUE, t_prime = t_prime)

survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0:strata(t_group) + A0:strata(t_group)
                  + N2:strata(t_group) + N3:strata(t_group), data = data)
survfit1 <- coxph(Surv(tstart, tstop, Delta == 1)~ L0:strata(t_group) + A0:strata(t_group)
                  + N2:strata(t_group) + N3:strata(t_group), data = data)
survfit2 <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0:strata(t_group) + A0:strata(t_group)
                  + N2:strata(t_group) + N3:strata(t_group), data = data)
survfit3 <- coxph(Surv(tstart, tstop, Delta == 3)~ L0:strata(t_group) + A0:strata(t_group)
                  + N2:strata(t_group) + N3:strata(t_group), data = data)


any(abs(survfit0$coefficients[c(1,3,5,7)] - beta[c(1,2,5,6),1]) > 0.1)
any(abs(survfit1$coefficients[c(1,3,5,7)] - beta[c(1,2,5,6),2]) > 0.1)
any(abs(survfit2$coefficients[c(1,3,5,7)] - beta[c(1,2,5,6),3]) > 0.1)
any(abs(survfit3$coefficients[c(1,3,5,7)] - beta[c(1,2,5,6),4]) > 0.1)

betap <- beta + time_var_eff
any(abs(survfit0$coefficients[c(2,4,6,8)] - betap[c(1,2,5,6),1]) > 0.1)
any(abs(survfit1$coefficients[c(2,4,6,8)] - betap[c(1,2,5,6),2]) > 0.1)
any(abs(survfit2$coefficients[c(2,4,6,8)] - betap[c(1,2,5,6),3]) > 0.1)
any(abs(survfit3$coefficients[c(2,4,6,8)] - betap[c(1,2,5,6),4]) > 0.1)

# Time varying effects and tests: simDisease

data <- simDisease(N = 10^4, eta = rep(0.1,3), nu = rep(1.1,3),  cens = 1,
                   beta_L0_D = 1, beta_L0_L = 1, beta_L_D = 1, beta_A0_D = 0,
                   beta_A0_L = 0, beta_L0_C = 0, beta_A0_C = 0, beta_L_C = 0,
                   followup = Inf, lower = 10^(-15), upper = 200,
                   beta_L_D_t_prime = 1, t_prime = t_prime, gen_A0 = NULL)

data <- IntFormatData(data, timeVar = TRUE, t_prime = t_prime, N_cols = 6)

survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0:strata(t_group) + A0:strata(t_group)
                  + L:strata(t_group), data = data)
survfit1 <- coxph(Surv(tstart, tstop, Delta == 1)~ L0:strata(t_group) + A0:strata(t_group)
                  + L:strata(t_group), data = data)
survfit2 <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0:strata(t_group) + A0:strata(t_group), 
                  data = data[L == 0])

# Time varying effects and tests: simDropIn

data <- simDropIn(N = 10^4, beta_L_A = 1, beta_L_Z = 2, beta_L_D = 1.5, beta_L_C = 0,
                  beta_A_L = -0.5,  beta_A_Z = -0.5, beta_A_D = -1, beta_A_C = 0,
                  beta_Z_L = -1, beta_Z_A = 0, beta_Z_D = -1, beta_Z_C = 0,
                  beta_L0_L = 1, beta_L0_A = 1, beta_L0_Z = 1, beta_L0_D = 1, beta_L0_C = 0,
                  beta_A0_L = -1.5, beta_A0_A = 0, beta_A0_Z = 0, beta_A0_D = -2, beta_A0_C = 0,
                  eta = c(0.5, 0.5, 0.1, 0.25),  nu = c(1.1, 1.1, 1.1, 1.1),
                  adherence = FALSE, followup = Inf, cens = 1,
                  generate.A0 = function(N, L0) stats::rbinom(N, 1, 0.5),
                  lower = 1e-200, upper = 1e10,
                  beta_L_A_prime = 1, beta_L_Z_prime = 0, beta_L_D_prime = 0,
                  beta_L_C_prime = 1, beta_A_L_prime = 0, beta_A_Z_prime = 0,
                  beta_A_D_prime = 1, beta_A_C_prime = 0, beta_Z_L_prime = 0,
                  beta_Z_A_prime = 1, beta_Z_D_prime = 0, beta_Z_C_prime = 0,
                  beta_L0_L_prime = 1,beta_L0_A_prime = 0,beta_L0_Z_prime = 0,
                  beta_L0_D_prime = 1,beta_L0_C_prime = 0,beta_A0_L_prime = 0,
                  beta_A0_A_prime = 0,beta_A0_Z_prime = 0,beta_A0_D_prime = 0,
                  beta_A0_C_prime = 0, t_prime = 1)

data <- IntFormatData(data, timeVar = TRUE, t_prime = t_prime, N_cols = 6:7)

survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0:strata(t_group) + A0:strata(t_group)
                  + Z:strata(t_group) + L:strata(t_group), data = data)
survfit1 <- coxph(Surv(tstart, tstop, Delta == 1)~ L0:strata(t_group) + A0:strata(t_group)
                  + Z:strata(t_group) + L:strata(t_group), data = data)
survfit2 <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0:strata(t_group) + A0:strata(t_group)
                  +  L:strata(t_group), data = data[L == 0])
survfit3 <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0:strata(t_group) + A0:strata(t_group)
                  +  Z:strata(t_group), data = data[L == 0])
