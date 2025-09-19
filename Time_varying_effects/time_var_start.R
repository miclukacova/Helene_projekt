# Time varying effects and tests

# Tjek af at det virker når vi ikke har time varying effects

eta <- rep(0.1, 4)
term_deltas <- c(0,1)
time_var_eff <- matrix(0.1, ncol = 4, nrow = 6)
beta <- matrix(nrow = 6, ncol = 4, rnorm(4*6, sd = 0.5))
# t_prime bliver sat højt så vi altid er i det sædvanlige scenarie
t_prime <- 1000                                                

data <- simEventDataTimeVar(N = 2*10^4, t_prime = t_prime, time_var_eff = time_var_eff, eta = eta,
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
data <- simEventDataTimeVar(N = 2*10^4, t_prime = t_prime, time_var_eff = time_var_eff, eta = eta,
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

data <- simEventDataTimeVar(N = 2*10^5, t_prime = t_prime, time_var_eff = time_var_eff, eta = eta,
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

