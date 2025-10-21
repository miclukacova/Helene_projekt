at_risk_cov <- function(covariates) {
  return(c(
    1,1,
    as.numeric(covariates[1] < 0.5),
    as.numeric(covariates[1] < 0.5)))
}

# Checking the various functions

#simEventData
data <- simEventData(N = 20000, at_risk_cov = at_risk_cov, beta = matrix(1, nrow = 6, ncol = 4))
data_int <- IntFormatData(data)
survfit0 <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A0 + N2 + N3, 
                  data = data_int[L0 < 0.5])

# simTreatment
data <- simTreatment(N = 20000, at_risk_cov = at_risk_cov)
data_int <- IntFormatData(data, N_cols = c(5,6))
survfit0 <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + A, 
                  data = data_int[(L0 < 0.5) & (L < 1)])

# simEventTV
data <- simEventTV(N = 20000, at_risk_cov = at_risk_cov, beta = matrix(1, nrow = 6, ncol = 4))
data_int <- IntFormatData(data)
survfit0 <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + N2 + N3, 
                  data = data_int[A0 > 0.5])

# simDropIn
data <- simDropIn(N = 20000, at_risk_cov = at_risk_cov)
data_int <- IntFormatData(data, N_cols = c(6,7))
survfit0 <- coxph(Surv(tstart, tstop, Delta == 3) ~ A0 + L0 + Z, 
                  data = data_int[L0 < 0.5 & L < 1])
# simDisease
at_risk_cov <- function(covariates) {
  return(c(
    1,1,
    as.numeric(covariates[1] < 0.5)))
}
data <- simDisease(N = 20000, at_risk_cov = at_risk_cov)
data_int <- IntFormatData(data, N_cols = c(6))
survfit0 <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A0, 
                  data = data_int[L0 < 0.5 & L < 1])


## Interaktionseffekter

set.seed(736)
beta <- matrix(0, ncol = 3, nrow = 5)

at_risk <- function(events) {
  return(c(
    1,1,                         # Always at risk for event 0 and 1
    as.numeric(events[3] < 2)))   # Event 2 can occur twice
}


data <- simEventData(N = 10000, beta = beta, at_risk = at_risk,
                     override_beta =list("N2 * A0" = c("N1" = 2)))

data_int <- IntFormatData(data, N_cols = 6:8)

survfit0 <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + N2+ N2:A0, 
                  data = data_int)
