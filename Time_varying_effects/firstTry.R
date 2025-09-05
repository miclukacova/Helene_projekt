###### Try time varying effects
library(simevent)

eta <- rep(0.1, 3)
term_deltas <- c(0,1,2)
time_var_eff <- matrix(1, ncol = 3, nrow = 5)

datatry <- simEventDataTry(N = 10000, t_prime = 1, time_var_eff = time_var_eff, eta = eta,
                term_deltas = term_deltas)

#plotEventData(datatry[Time <= 1])
#plotEventData(datatry[Time > 1])

datatry[, A0_prime := (Time > 1) * A0]
datatry[, L0_prime := (Time > 1) * L0]

survfit0 <- coxph(Surv(Time, Delta == 0) ~ L0 + A0 + L0_prime + A0_prime, data = datatry)
survfit1 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0 + L0_prime + A0_prime, data = datatry)
survfit2 <- coxph(Surv(Time, Delta == 2) ~ L0 + A0 + L0_prime + A0_prime, data = datatry)
