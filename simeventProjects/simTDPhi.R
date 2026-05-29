################### Trying out simEventDataTdPhi ###############################

# To do's:
# * store værdier af beta2

beta <- matrix(rnorm(4*6, 0, 0.1), ncol = 4)

data <- simEventDataTdPhi(N = 10^4, beta2 = c(0.01,2,0,-0.1),
                          beta = beta, upper = 10^20, max_iter = 10^3)

# Sanity check
#data <- simEventDataTdPhi(N = 1000, beta2 = c(0,0,0,0),
#                          beta = beta, upper = 10^10, max_iter = 500)

# Check Simulations - do they work?

# Transform data into tstart tstop format
data_int <- IntFormatData(data)

# Fit models
survfit_cens <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + N2 + N3, data = data_int)
survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + N2 + N3, data = data_int)
survfit_oper <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A0 + N2 + N3, data = data_int)
survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + A0 + N2 + N3, data = data_int)

# Compare confidence intervals and true values
expect_true(confint(survfit_cens, level = 0.99)[1,1] <= beta[1,1] & beta[1,1] <= confint(survfit_cens, level = 0.95)[1,2])
expect_true(confint(survfit_cens, level = 0.99)[2,1] <= beta[2,1] & beta[2,1] <= confint(survfit_cens, level = 0.95)[2,2])
expect_true(confint(survfit_cens, level = 0.95)[3,1] <= beta[5,1] & beta[5,1] <= confint(survfit_cens, level = 0.95)[3,2])
expect_true(confint(survfit_cens, level = 0.95)[4,1] <= beta[6,1] & beta[6,1] <= confint(survfit_cens, level = 0.95)[4,2])
expect_true(confint(survfit_death, level = 0.95)[1,1] <= beta[1,2] & beta[1,2] <= confint(survfit_death, level = 0.95)[1,2])
expect_true(confint(survfit_death, level = 0.95)[2,1] <= beta[2,2] & beta[2,2] <= confint(survfit_death, level = 0.95)[2,2])
expect_true(confint(survfit_death, level = 0.95)[3,1] <= beta[5,2] & beta[5,2] <= confint(survfit_death, level = 0.95)[3,2])
expect_true(confint(survfit_death, level = 0.95)[4,1] <= beta[6,2] & beta[6,2] <= confint(survfit_death, level = 0.95)[4,2])
expect_true(confint(survfit_oper, level = 0.95)[1,1] <= beta[1,3] & beta[1,3] <= confint(survfit_oper, level = 0.95)[1,2])
expect_true(confint(survfit_oper, level = 0.95)[2,1] <= beta[2,3] & beta[2,3] <= confint(survfit_oper, level = 0.95)[2,2])
expect_true(confint(survfit_oper, level = 0.95)[3,1] <= beta[5,3] & beta[5,3] <= confint(survfit_oper, level = 0.95)[3,2])
expect_true(confint(survfit_oper, level = 0.95)[4,1] <= beta[6,3] & beta[6,3] <= confint(survfit_oper, level = 0.95)[4,2])
expect_true(confint(survfit_cov, level = 0.95)[1,1] <= beta[1,4] & beta[1,4] <= confint(survfit_cov, level = 0.95)[1,2])
expect_true(confint(survfit_cov, level = 0.95)[2,1] <= beta[2,4] & beta[2,4] <= confint(survfit_cov, level = 0.95)[2,2])
expect_true(confint(survfit_cov, level = 0.95)[3,1] <= beta[5,4] & beta[5,4] <= confint(survfit_cov, level = 0.95)[3,2])
expect_true(confint(survfit_cov, level = 0.95)[4,1] <= beta[6,4] & beta[6,4] <= confint(survfit_cov, level = 0.95)[4,2])

# Jeg tror det virker????
