at_risk <- function(events, covariates) {
  return(c(
    1,1,
    as.numeric(events[3] < 2 & (covariates[1] < 0.5)),
    as.numeric(events[4] < 1 & (covariates[2] > 0.5))))
}

at_risk(c(1,1,1,1), c(0.1,0.1))

sim_data <- simEventData(N = 10, at_risk = at_risk)
