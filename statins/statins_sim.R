#-------------------------------------------------------------------------------
# Description
#-------------------------------------------------------------------------------

# This script first simulates a data set from the statin setting
# Then beta, coefficients are estimated and compared to the inputted
# Then the function simeventcox is used to simulate new data from the fitted cox models. The estimate coefficients 
# are compared to the original coefficients.
# The last part of the script performs interventions and calculates interventions effects. 

#-------------------------------------------------------------------------------
# Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(simevent)
library(survival)
library(ggplot2)
library(ggpubr)

#-------------------------------------------------------------------------------
# Simulate trial data set
#-------------------------------------------------------------------------------
# There are 12 processes: Censoring, Death, CVD, Off Statins, On Statins, Treatment, Disease, LDL 0, LDL 1, LDL 2, LDL 3, LDL 4

n_cov <- 6
n_proc <- 12

# Regression coefficients

# Creating beta matrix
beta <- matrix(0, nrow = n_cov + n_proc, ncol = n_proc)
rownames(beta) <- c(
  "L0", "A0",
  "L1", "L2", "L3", "L4", 
  "C", "D", "CVD", "S0", "S1", "L", "A", "LDL0", "LDL1", "LDL2", "LDL3", "LDL4"
)

colnames(beta) <- c("C", "D", "CVD", "S0", "S1", "L", "A", "LDL0", "LDL1", "LDL2", "LDL3", "LDL4")


# The coefficients cannot be specified without the simulaitons breaking down.  
#beta[sample(nrow(beta)*ncol(beta), 30)] <- rep(seq(-0.3,0.3, length.out = 15),2)


#beta[c(1:n_cov,(n_cov + 4):nrow(beta)),2] <- c(-0.442, 0.065, -0.157, -0.571, -0.131, -0.242, -0.112,
#                                      0.077, 0.371, 0.559, 0.479, -0.211, 0.140)
#
#beta[c(1:n_cov,(n_cov + 4):nrow(beta)),3] <- c(-0.592, 0.035, 0.774, 0.182, 0.555,
#                                      -0.027, 0.029, 0.040, 0.071, 0.977,
#                                       0.153, 0.180, 0.148)
#
#beta[c(1:n_cov,(n_cov + 4):nrow(beta)),4] <- c(0.013, -0.052, 0.029, -0.009, 0.073,
#                                     -0.028, 0.039, 0.031, 0, 0.056,
#                                      0.031, 0.028, 0.004)
#
#beta[c(1:n_cov,(n_cov + 4):nrow(beta)),5] <- c(-0.253, -0.016, 0, 0.005, 0.143,
#                                       0.450, -0.015, 0.042, -0.252,
#                                       0.036, 0.281, 0.079, 0.208)
#
#beta[c(1:n_cov,(n_cov + 4):nrow(beta)),6] <- c(-0.031, -0.058, 0.092, 0.053, 0.199,
#                                      -0.036, -0.001, 0.061, 0.006, 0.339,
#                                       0.024, 0.057, 0.129)
#
#beta[c(1:n_cov,(n_cov + 4):nrow(beta)),7] <- c(0.128, -0.074, 0.113, 0.190, 0.289,
#                                     -0.047, -0.384, -0.015, -0.162, 0.162,
#                                     -0.021, -0.331, 0.619)
#
#beta[c(1:n_cov,(n_cov + 4):nrow(beta)),8] <- c(-0.183, 0.061, 0.061, 0.063, 0.214,
#                                      0.021, 0.106, 0.043, -0.503,
#                                      0.142, 0.041, 0.775, -0.361)
# Covariate generating distribution
add_cov <- list()

gen_L0 <- function(N) rbinom(N, 1, 0.4)                                         # koen
gen_A0 <- function(N, L0) pmin(rexp(N, 0.3) + 70, 100)                          # alder
add_cov[[1]] <- function(N) sample(c(1,2,3,4), N, replace = TRUE)               # base_civst
add_cov[[2]] <- function(N) rpois(N, 0.25)                                      # n_diag_base
add_cov[[3]] <- function(N) pmax(rnorm(N, 2, 1), 0.2)                           # base_LDL
add_cov[[4]] <- function(N) rpois(N, 5)                                         # base_drugs

# at risk function
at_risk <- function(events) {
  return(c(1,                                                                   # Censoring
           1,                                                                   # If you have not you died yet are at risk
           1,                                                                   # If you have not had CVD you are at risk
           as.numeric((events[4] == events[5]) & (events[4] < 10)),             # You are at risk of Statin Stop if you have stopped as many times as you've started
           as.numeric((events[4] > events[5]) & (events[5] < 10)),              # You are at risk of Statin Start if you have stopped more times than you've started
           as.numeric(events[6] <= 10),                                         # Increase in number of diseases
           as.numeric(events[7] <= 10),                                         # Increase in number of medicines
           as.numeric(events[8] <= 10),                                         # LDL jump  (complicated at risk structure modelled later)
           as.numeric(events[9] <= 10),                                         # LDL jump  (complicated at risk structure modelled later)
           as.numeric(events[10] <= 10),                                        # LDL jump  (complicated at risk structure modelled later)
           as.numeric(events[11] <= 10),                                        # LDL jump  (complicated at risk structure modelled later)
           as.numeric(events[12] <= 10)))                                       # LDL jump  (complicated at risk structure modelled later)
  
}


# Estimerede parametre, som vi har fået ved at fitte lm fits
#nu <- c(1, 0.8528, 0.7549, 1.2640, 0.7433, 0.8939, 0.7158, 0.6953)
#eta <- c(1, 0.0029, 0.0005, 0.0065, 0.0229, 0.0561, 0.0432, 0.0497)

# Estimerede parametre, som vi har fået ved at fitte på event of interest og terminale events
#eta <- c(1, 0.0007, 0.0005, 0.0055, 0.0170, 0.0462, 0.0674, 0.0477)
#nu <- c(1, 1.2383, 0.8755, 1.3726, 0.7572, 0.8640, 0.6226, 0.6768)

# temporary adjustments to get reasonable simulations
#eta[c(1,2,8)] <- eta[c(1,2,8)] / 3
#nu[c(6,7)] <- nu[c(6,7)] * 2

eta <- rep(0.001, n_proc)
nu <- rep(1.001, n_proc)
eta[5] <- 0.1

# Simulating from simStatinData
data <- simStatinData(beta = beta, 
                      N = 10^5, 
                      add_cov = add_cov, 
                      followup = 60, 
                      gen_A0 = gen_A0, 
                      gen_L0 = gen_L0,
                      eta = eta, 
                      nu = nu,
                      lower = 10^(-20),
                      upper = 10^4,
                      at_risk = at_risk)

# We name the processes
colnames(data)[(6 + length(add_cov)):ncol(data)] <- colnames(beta) 

plotEventData(data[1:2000,])
data <- IntFormatData(data, N_cols = (n_cov + 4):(n_cov+n_proc+3))

#-------------------------------------------------------------------------------
# Estimating beta coefficiens
#-------------------------------------------------------------------------------

# Fit models
# Models where the indidividuals are always at risk
vars <- c("L0", "A0", "L1", "L2", "L3", "L4", colnames(beta)[4:n_proc])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 0) ~", paste(vars, collapse = " + ")))
survfit0 <- coxph(form, data = data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 1) ~", paste(vars, collapse = " + ")))
survfit1 <- coxph(form, data = data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 2) ~", paste(vars, collapse = " + ")))
survfit2 <- coxph(form, data = data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 3) ~", paste(vars, collapse = " + ")))
survfit3 <- coxph(form, data = data[S0 == S1])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 4) ~", paste(vars, collapse = " + ")))
survfit4 <- coxph(form, data = data[S0 > S1])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 5) ~", paste(vars, collapse = " + ")))
survfit5 <- coxph(form, data = data[L <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 6) ~", paste(vars, collapse = " + ")))
survfit6 <- coxph(form, data = data[A <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 7) ~", paste(vars, collapse = " + ")))
survfit7 <- coxph(form, data = data[LDL0 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 8) ~", paste(vars, collapse = " + ")))
survfit8 <- coxph(form, data = data[LDL1 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 9) ~", paste(vars, collapse = " + ")))
survfit9 <- coxph(form, data = data[LDL2 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 10) ~", paste(vars, collapse = " + ")))
survfit10 <- coxph(form, data = data[LDL3 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 11) ~", paste(vars, collapse = " + ")))
survfit11 <- coxph(form, data = data[LDL4 <= 10])

# Check estimations
beta_comp <- beta[c(1:n_cov,((n_cov + 4):nrow(beta))), ]

df <- rbind(cbind(beta_comp[,1],  confint(survfit0), 0), cbind(beta_comp[,2],  confint(survfit1), 1))
df <- rbind(df, cbind(beta_comp[,3],  confint(survfit2), 2))
df <- rbind(df, cbind(beta_comp[,4],  confint(survfit3), 3))
df <- rbind(df, cbind(beta_comp[,5],  confint(survfit4), 4))
df <- rbind(df, cbind(beta_comp[,6],  confint(survfit5), 5))
df <- rbind(df, cbind(beta_comp[,7],  confint(survfit6), 6))
df <- rbind(df, cbind(beta_comp[,8],  confint(survfit7), 7))
df <- rbind(df, cbind(beta_comp[,9],  confint(survfit8), 8))
df <- rbind(df, cbind(beta_comp[,10],  confint(survfit9), 9))
df <- rbind(df, cbind(beta_comp[,11],  confint(survfit5), 10))
df <- rbind(df, cbind(beta_comp[,12],  confint(survfit5), 11))

colnames(df) <- c("actual", "LowCI", "UpCI", "fit")
df <- data.table(df)
df[, Index := seq_len(.N), by = fit]
df[, In := actual >= LowCI & actual <= UpCI]


ggplot(df, aes(x = actual, y = Index)) +
  geom_segment(aes(x = LowCI, xend = UpCI, yend = Index),
               colour = "grey50", linewidth = 2) +
  geom_point(aes(color = In), size = 2.5) +
  facet_wrap(~fit, scales = "free") +
  theme_bw()+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "darkgreen"))

#-------------------------------------------------------------------------------
# Simulating new data without intervention
#-------------------------------------------------------------------------------

old_vars <-  data[,4:(3+n_cov)]

cox_fits <- list(
  C = survfit0,
  D = survfit1,
  CVD = survfit2, 
  S0 = survfit3,
  S1 =survfit4, 
  L = survfit5, 
  A = survfit6,
  LDL0 = survfit7, 
  LDL1 = survfit8, 
  LDL2 = survfit9, 
  LDL3 = survfit10, 
  LDL4 = survfit11 
)


sim_data0 <- simEventCox(
  N = 10^4,
  cox_fits = cox_fits,
  old_vars = old_vars,
  n_event_max = c(1, 1, 1, rep(10, 9)),
  term_events = c(1, 2, 3),
  at_risk = at_risk
)

xlimm <- max(range(data[1:500,Time]),range(sim_data0[1:500,Time]))
ggarrange(plotEventData(sim_data0[1:500,], title = "Simulated") + xlim(c(0,xlimm)), 
          plotEventData(data[1:500,], title = "Original Data")+ xlim(c(0,xlimm)), ncol = 2)


# Sanity check of whether the simulated data corresponds to the original data
sim_data <- IntFormatData(sim_data0, N_cols = (n_cov + 4):(n_cov+n_proc+3))

# Fit models
# Models where the indidividuals are always at risk
vars <- setdiff(names(sim_data), c("ID", "Time", "k", "C", "D", "CVD", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 0) ~", paste(vars, collapse = " + ")))
survfit0 <- coxph(form, data = sim_data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 1) ~", paste(vars, collapse = " + ")))
survfit1 <- coxph(form, data = sim_data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 2) ~", paste(vars, collapse = " + ")))
survfit2 <- coxph(form, data = sim_data)

# Models for which the individual is not always at risk
form <- as.formula(paste("Surv(tstart, tstop, Delta == 3) ~", paste(vars, collapse = " + ")))
survfit3 <- coxph(form, data = sim_data[S0 == S1])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 4) ~", paste(vars, collapse = " + ")))
survfit4 <- coxph(form, data = sim_data[S0 > S1])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 5) ~", paste(vars, collapse = " + ")))
survfit5 <- coxph(form, data = sim_data[L <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 6) ~", paste(vars, collapse = " + ")))
survfit6 <- coxph(form, data = sim_data[A <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 7) ~", paste(vars, collapse = " + ")))
survfit7 <- coxph(form, data = sim_data[LDL0 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 8) ~", paste(vars, collapse = " + ")))
survfit8 <- coxph(form, data = sim_data[LDL1 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 9) ~", paste(vars, collapse = " + ")))
survfit9 <- coxph(form, data = sim_data[LDL2 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 10) ~", paste(vars, collapse = " + ")))
survfit10 <- coxph(form, data = sim_data[LDL3 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 11) ~", paste(vars, collapse = " + ")))
survfit11 <- coxph(form, data = sim_data[LDL4 <= 10])


# Check estimations
beta_comp <- beta[c(1:n_cov,((n_cov + 4):nrow(beta))), ]

df <- rbind(cbind(beta_comp[,1],  confint(survfit0), 0), cbind(beta_comp[,2],  confint(survfit1), 1))
df <- rbind(df, cbind(beta_comp[,3],  confint(survfit2), 2))
df <- rbind(df, cbind(beta_comp[,4],  confint(survfit3), 3))
df <- rbind(df, cbind(beta_comp[,5],  confint(survfit4), 4))
df <- rbind(df, cbind(beta_comp[,6],  confint(survfit5), 5))
df <- rbind(df, cbind(beta_comp[,7],  confint(survfit6), 6))
df <- rbind(df, cbind(beta_comp[,8],  confint(survfit7), 7))
df <- rbind(df, cbind(beta_comp[,9],  confint(survfit8), 8))
df <- rbind(df, cbind(beta_comp[,10],  confint(survfit9), 9))
df <- rbind(df, cbind(beta_comp[,11],  confint(survfit10), 10))
df <- rbind(df, cbind(beta_comp[,12],  confint(survfit11), 11))

colnames(df) <- c("actual", "LowCI", "UpCI", "fit")
df <- data.table(df)
df[, Index := seq_len(.N), by = fit]
df[, In := actual >= LowCI & actual <= UpCI]

ggplot(df, aes(x = actual, y = Index)) +
  geom_segment(aes(x = LowCI, xend = UpCI, yend = Index),
               colour = "grey50", linewidth = 2) +
  geom_point(aes(color = In), size = 2.5) +
  facet_wrap(~fit, scales = "free") +
  theme_bw()+
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "darkgreen"))


#-------------------------------------------------------------------------------
# Simulating new data with intervention
#-------------------------------------------------------------------------------

# An intervention multiplying the off statins process by alpha
alpha <- 2
intervention <- function(j, basehaz) if(j == 4) alpha * basehaz else basehaz

sim_data_int <- simEventCox(
  N = 10^4,
  cox_fits = cox_fits,
  old_vars = old_vars,
  n_event_max = c(1, 1, 1, rep(10, 9)),
  term_events = c(1, 2, 3),
  at_risk = at_risk,
  intervention2 = intervention
)

ggarrange(plotEventData(sim_data0[ID %in% 1:200], title = "Non intervened"), 
          plotEventData(sim_data_int[ID %in% 1:200], title = "Intervened"), ncol = 2)

#-------------------------------------------------------------------------------
# Calculate Intervention Effects
#-------------------------------------------------------------------------------

tau <- 5

# Proportion of subjects dying before time $\tau$
mean(sim_data0[, any(Delta == 1 & Time < tau)[1], by = "ID"][[2]])
mean(sim_data_int[, any(Delta == 1 & Time < tau)[1], by = "ID"][[2]])

# Proportion of subjects experiencing MACE before time $\tau$
mean(sim_data0[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]])
mean(sim_data_int[, any(Delta == 2 & Time < tau)[1], by = "ID"][[2]])

# Proportion of subjects experiencing Disease before time tau
mean(sim_data0[, any(Delta == 5 & Time < tau)[1], by = "ID"][[2]])
mean(sim_data_int[, any(Delta == 5 & Time < tau)[1], by = "ID"][[2]])

# Proportion of subjects experiencing Treatment before time tau
mean(sim_data0[, any(Delta == 4 & Time < tau)[1], by = "ID"][[2]])
mean(sim_data_int[, any(Delta == 4 & Time < tau)[1], by = "ID"][[2]])

#-------------------------------------------------------------------------------
# Plot of "5-year risk of MACE" and "5-year risk of death"
#-------------------------------------------------------------------------------

# at risk function is modified so that you are no longer at risk for censoring
at_risk <- function(events) {
  return(c(0,                                                                   # Censoring
           1,                                                                   # If you have not you died yet are at risk
           1,                                                                   # If you have not had CVD you are at risk
           as.numeric((events[4] == events[5]) & (events[4] < 10)),             # You are at risk of Statin Stop if you have stopped as many times as you've started
           as.numeric((events[4] > events[5]) & (events[5] < 10)),              # You are at risk of Statin Start if you have stopped more times than you've started
           as.numeric(events[6] <= 10),                                         # Increase in number of diseases
           as.numeric(events[7] <= 10),                                         # Increase in number of medicines
           as.numeric(events[8] <= 10),                                         # LDL jump  (complicated at risk structure modelled later)
           as.numeric(events[9] <= 10),                                         # LDL jump  (complicated at risk structure modelled later)
           as.numeric(events[10] <= 10),                                        # LDL jump  (complicated at risk structure modelled later)
           as.numeric(events[11] <= 10),                                        # LDL jump  (complicated at risk structure modelled later)
           as.numeric(events[12] <= 10)))                                       # LDL jump  (complicated at risk structure modelled later)
  
}

# The vector of alphas for the alpha interventions
alphas <- seq(0.5,2, by = 0.1)
# The r
risk_alpha <- matrix(nrow = length(alphas), ncol = 2)
  
for(i in seq_along(alphas)){
  print(i)
  res_sim <- alphaSim(N = 1e5,
                      eta =  eta,
                      nu = nu,
                      alpha = alphas[i],
                      tau = 5,
                      setting = "Statin",
                      beta = beta, 
                      add_cov = add_cov)
   
   risk_alpha[i,] <- unlist(res_sim)
}

risk_alpha <- data.frame(risk_alpha)
colnames(risk_alpha) <- c("Death", "MACE")

ggplot(risk_alpha, aes(x = alphas, y = Death)) +
  geom_line(aes(color = "Death"), linewidth = 2) +
  geom_line(aes(y = MACE, color = "MACE"), linewidth = 2) +
  theme_bw()+
  scale_color_manual(values = c("MACE" = "red", "Death" = "darkgreen"))+
  xlab(expression(alpha))+
  ylab("Risk")


