#-------------------------------------------------------------------------------
# Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(simevent)
library(survival)
library(ggplot2)
library(ggpubr)

#-------------------------------------------------------------------------------
# Simulating a dummy trial data set
#-------------------------------------------------------------------------------

# Regression coefficients
#beta <- matrix(0, nrow = 30, ncol = 6)
beta <- matrix(runif(30*6, -1, 1), nrow = 30, ncol = 6)

# Covariate generating distribution
data0 <- matrix(sample(c(runif(100*22), rep(0,100*22)), 200*22), nrow = 200, ncol = 22)
add_cov <- list()
for(i in 1:22){
  add_cov[[i]] <- function(N) sample(data0[,i], N, replace = TRUE)
}

# Simulating from simStatinData
data <- simStatinData(beta = beta, N = 3*10^4, add_cov = add_cov, followup = Inf)
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

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 4) ~", paste(vars, collapse = " + ")))
survfit5 <- coxph(form, data = data[A <= 3])

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD",  "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 5) ~", paste(vars, collapse = " + ")))
survfit6 <- coxph(form, data = data[L <= 3])

# Check estimations
df <- rbind(cbind(beta[c(1:24,28:30),1],  confint(survfit1), 1),
                  cbind(beta[c(1:24,28:30),2],  confint(survfit2), 2))

df <- rbind(df, cbind(beta[c(1:24,28:30),3],  confint(survfit3), 3))
df <- rbind(df, cbind(beta[c(1:24,29:30),4],  confint(survfit4), 4))
df <- rbind(df, cbind(beta[c(1:24,28:30),5],  confint(survfit5), 5))
df <- rbind(df, cbind(beta[c(1:24,28:30),6],  confint(survfit6), 6))



colnames(df) <- c("actual", "LowCI", "UpCI", "fit")
df <- data.table(df)
df[, Index := seq_len(.N), by = fit]
df[, In := actual >= LowCI & actual <= UpCI]

ggplot(df, aes(x = actual, y = Index)) +
  geom_segment(aes(x = LowCI, xend = UpCI, yend = Index),
               colour = "grey50", linewidth = 2) +
  geom_point(aes(color = In), size = 2.5) +
  facet_wrap(~fit) +
  theme_bw()+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "darkgreen"))

#-------------------------------------------------------------------------------
# Simulating new data without intervention
#-------------------------------------------------------------------------------

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

sim_data0 <- simEventCox(
  10^4,
  cox_fits,
  list_old_vars = list_old_vars,
  n_event_max = c(1, 1, 1, 1, 3, 3),
  term_events = c(1, 2, 3),
)


ggarrange(plotEventData(sim_data0[1:500,], title = "Simulated"), 
          plotEventData(data[1:500,], title = "Original Data"), ncol = 2)


# Sanity check of whether the simulated data corresponds to the original data
sim_data <- IntFormatData(sim_data0, N_cols = 28:33)
#sim_data[,Delta := Delta -1]

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

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 4) ~", paste(vars, collapse = " + ")))
survfit5 <- coxph(form, data = sim_data[A < 3])

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 5) ~", paste(vars, collapse = " + ")))
survfit6 <- coxph(form, data = sim_data[L < 3])

# Check estimations
df <- rbind(cbind(beta[c(1:24,28:30),1],  confint(survfit1), 1),
            cbind(beta[c(1:24,28:30),2],  confint(survfit2), 2))

df <- rbind(df, cbind(beta[c(1:24,28:30),3],  confint(survfit3), 3))
df <- rbind(df, cbind(beta[c(1:24,29:30),4],  confint(survfit4), 4))
df <- rbind(df, cbind(beta[c(1:24,28:30),5],  confint(survfit5), 5))
df <- rbind(df, cbind(beta[c(1:24,28:30),6],  confint(survfit6), 6))



colnames(df) <- c("actual", "LowCI", "UpCI", "fit")
df <- data.table(df)
df[, Index := seq_len(.N), by = fit]
df[, In := actual >= LowCI & actual <= UpCI]

ggplot(df, aes(x = actual, y = Index)) +
  geom_segment(aes(x = LowCI, xend = UpCI, yend = Index),
               colour = "grey50", linewidth = 2) +
  geom_point(aes(color = In), size = 2.5) +
  facet_wrap(~fit) +
  theme_bw()+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "darkgreen"))


#-------------------------------------------------------------------------------
# Simulating new data with intervention
#-------------------------------------------------------------------------------

# An intervention multiplying the treatment process by alpha
alpha <- 2
intervention <- function(j, basehaz) if(j == 5) alpha * basehaz else basehaz

sim_data_int <- simEventCox(
  10^4,
  cox_fits,
  list_old_vars = list_old_vars,
  n_event_max = c(1, 1, 1, 1, 3, 3),
  term_events = c(1, 2, 3),
  intervention2 = intervention
)


ggarrange(plotEventData(sim_data0[ID %in% 1:200], title = "Non intervened"), 
          plotEventData(sim_data_int[ID %in% 1:200], title = "Intervened"), ncol = 2)



#-------------------------------------------------------------------------------
# Calculate Intervention Effects
#-------------------------------------------------------------------------------

tau <- 2

# Proportion of subjects dying before some time $\tau$
sim_data0[Delta == 1, mean(Delta == 1 & Time < tau)]
sim_data_int[Delta == 1, mean(Delta == 1 & Time < tau)]

# Proportion of subjects experiencing Disease
mean(sim_data0[, any(Delta == 5 & Time < tau)[1], by = "ID"][[2]])
mean(sim_data_int[, any(Delta == 5 & Time < tau)[1], by = "ID"][[2]])

# Proportion of subjects experiencing Treatment
mean(sim_data0[, any(Delta == 4 & Time < tau)[1], by = "ID"][[2]])
mean(sim_data_int[, any(Delta == 4 & Time < tau)[1], by = "ID"][[2]])


