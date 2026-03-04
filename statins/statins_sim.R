#-------------------------------------------------------------------------------
# Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(simevent)
library(survival)

#-------------------------------------------------------------------------------
# Simulating a dummy trial data set
#-------------------------------------------------------------------------------

# Regression coefficients
beta <- matrix(runif(30*6, -1, 1), nrow = 30, ncol = 6)

# Covariate generating distribution
data0 <- matrix(sample(c(runif(100*22), rep(0,100*22)), 200*22), nrow = 200, ncol = 22)
add_cov <- list()
for(i in 1:22){
  add_cov[[i]] <- function(N) sample(data0[,i], N, replace = TRUE)
}

# Simulating from simStatinData
data <- simStatinData(beta = beta, N = 10^6, add_cov = add_cov)
data <- IntFormatData(data, N_cols = 28:33)

#-------------------------------------------------------------------------------
# Fitting Models
#-------------------------------------------------------------------------------


# Fit models
# Models where the indidividuals are always at risk
vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 0) ~", paste(vars, collapse = " + ")))
survfit1 <- coxph(form, data = data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 1) ~", paste(vars, collapse = " + ")))
survfit2 <- coxph(form, data = data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 2) ~", paste(vars, collapse = " + ")))
survfit3 <- coxph(form, data = data)


# Models for which the individual is not always at risk
vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "OS", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 3) ~", paste(vars, collapse = " + ")))
survfit4 <- coxph(form, data = data[OS == 0])

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "A", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 4) ~", paste(vars, collapse = " + ")))
survfit5 <- coxph(form, data = data[A <= 3])

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "L", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 5) ~", paste(vars, collapse = " + ")))
survfit6 <- coxph(form, data = data[L <= 3])

# Check estimations
df <- rbind(rbind(cbind(beta[c(1:24,28:30),3],  confint(survfit3)),
                  cbind(beta[c(1:24,28:30),2],  confint(survfit2))), 
            cbind(beta[c(1:24,29:30),4],  confint(survfit4)))
colnames(df) <- c("actual", "LowCI", "UpCI")
df <- cbind(df, Index = 1:nrow(df))

ggplot(df) +
  geom_segment(aes(x = LowCI, xend = UpCI, 
                   y = Index, yend = Index),
               colour = "red") +
  geom_point(aes(x = actual, y = Index))


# Der går måske noget galt i process 0...???

#-------------------------------------------------------------------------------
# Simulating new data, respectively under intervention and without intervention
#-------------------------------------------------------------------------------

# Der er moget som går galt her....

list_old_vars <- list()
for(var in 4:27){
  list_old_vars[(var - 3)] <- data[,..var]
}

names(list_old_vars) <- colnames(data[, 4:27])

cox_fits <- list(
  "C" = survfit1,
  "D" = survfit2,
  "CVD" = survfit3,
  "OS" = survfit4,
  "A" = survfit5,
  "L" = survfit6
)


sim_data <- simEventCox(
  10^4,
  cox_fits,
  list_old_vars = list_old_vars,
  n_event_max = c(1, 1, 1, 1, 3, 3),
  term_events = c(1, 2, 3),
)

# Sanity check of whether the simulated data corresponds to the original data
sim_data <- IntFormatData(sim_data, N_cols = 28:33)

# Fit models
# Models where the indidividuals are always at risk
vars <- setdiff(names(sim_data), c("ID", "Time", "k", "C", "D", "CVD", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 0) ~", paste(vars, collapse = " + ")))
survfit1 <- coxph(form, data = sim_data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 1) ~", paste(vars, collapse = " + ")))
survfit2 <- coxph(form, data = sim_data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 2) ~", paste(vars, collapse = " + ")))
survfit3 <- coxph(form, data = sim_data)

# Models for which the individual is not always at risk
vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "OS", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 3) ~", paste(vars, collapse = " + ")))
survfit4 <- coxph(form, data = sim_data[OS == 0])

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "A", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 4) ~", paste(vars, collapse = " + ")))
survfit5 <- coxph(form, data = sim_data[A < 3])

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "L", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 5) ~", paste(vars, collapse = " + ")))
survfit6 <- coxph(form, data = sim_data[L < 3])

# Check estimations
df <- rbind(rbind(cbind(beta[c(1:24,28:30),1],  confint(survfit3)),
                  cbind(beta[c(1:24,28:30),2],  confint(survfit2))), 
            cbind(beta[c(1:24,29:30),4],  confint(survfit4)))
colnames(df) <- c("actual", "LowCI", "UpCI")
df <- cbind(df, Index = 1:nrow(df))

ggplot(df) +
  geom_segment(aes(x = LowCI, xend = UpCI, 
                   y = Index, yend = Index),
               colour = "red") +
  geom_point(aes(x = actual, y = Index))

