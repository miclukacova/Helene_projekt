data_tjek <- sim_event_data(N = 5 * 10^4, max_cens = 5)
data_tjek0 <- sim_event_data(N = 5 * 10^4)

# Transform data into tstart tstop format
data_int <- trans_int_data(data_tjek)
data_int0 <- trans_int_data(data_tjek0)

data_int[, at_risk_oper := as.numeric(A == 0)]
data_int[, at_risk_cov := as.numeric(L == 0)]
data_int0[, at_risk_oper := as.numeric(A == 0)]
data_int0[, at_risk_cov := as.numeric(L == 0)]

# Fit models
survfit_cens <- coxph(Surv(tstart, tstop, Delta == 0) ~ I(L0 / 50) + A + A0 + L,
                      data = data_int,
                      timefix = FALSE,
                      cluster = ID)
survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ I(L0 / 50) + A + A0 + L,
                       data = data_int,
                       timefix = FALSE,
                       cluster = ID)
survfit_oper <- coxph(Surv(tstart, tstop, Delta == 2) ~ I(L0 / 50) + A0 + L,
                      data = data_int[at_risk_oper == 1],
                      timefix = FALSE,
                      cluster = ID)
survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ I(L0 / 50) + A + A0,
                     timefix = FALSE,
                     data = data_int[at_risk_cov == 1],
                     cluster = ID)

# De her skulle gerne ikke have bias
survfit_cens0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ I(L0 / 50) + A + A0 + L,
                       data = data_int0,
                       timefix = FALSE,
                       cluster = ID)
survfit_death0 <- coxph(Surv(tstart, tstop, Delta == 1) ~ I(L0 / 50) + A + A0 + L,
                       data = data_int0,
                       timefix = FALSE,
                       cluster = ID)
survfit_oper0 <- coxph(Surv(tstart, tstop, Delta == 2) ~ I(L0 / 50) + A0 + L,
                       data = data_int0[at_risk_oper == 1],
                       timefix = FALSE,
                       cluster = ID)
survfit_cov0 <- coxph(Surv(tstart, tstop, Delta == 3) ~ I(L0 / 50) + A + A0,
                     timefix = FALSE,
                     data = data_int0[at_risk_cov == 1],
                     cluster = ID)

