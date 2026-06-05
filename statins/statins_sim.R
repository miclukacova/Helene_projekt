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
# There are 8 processes: Censoring, Death, CVD, Off Statins, Treatment, Disease, LDL increase, LDL decrease

n_cov <- 10
n_proc <- 8

# Regression coefficients

# Creating beta matrix
beta <- matrix(0, nrow = n_cov + n_proc, ncol = n_proc)
rownames(beta) <- c(
  "L0", "A0",
  "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8",
  "C", "D", "CVD", "OS", "L", "A", "LDL1", "LDL2"
)
colnames(beta) <- c("C", "D", "CVD", "OS", "L", "A", "LDL1", "LDL2")

#beta[c(1:10,14:nrow(beta)),1] <- c(-0.025, 0.001, 0.021, 0.005, -0.008, 
#                              -0.081, -0.008, 0.009, -0.023, 0.014,
#                              -0.047, -0.042, -0.024, -0.001, -0.058)

beta[c(1:10,14:nrow(beta)),2] <- c(-0.439, 0.065, 1.003, 0.599, 1.033,
                              1.198, -0.236, -0.110, 0, 0.077,
                              0.245, 0.557, 0.479, -0.201, 0.139)

# hvis jeg halvere koefficienterne her så ser det lidt bedre ud(?)
beta[,2] <- beta[,2]/4

beta[c(1:10,14:nrow(beta)),3] <- c(-0.558, 0.034, 0.776, 0.188, 0.558,
                              0, -0.024, 0.030, 0.129, 0.038,
                              -0.023, 0.974, 0.15, 0.186, 0.149)

beta[c(1:10,14:nrow(beta)),4] <- c(-0.007, -0.047, -0.007, -0.058, 0.040,
                              -0.032, -0.044, 0.030, 0, 0.033,
                              0, 0.099, 0.047, 0.011, 0.008)

beta[c(1:10,14:nrow(beta)),5] <- c(-0.251, -0.017, -0.423, -0.418, -0.278,
                              -0.449, 0.452, -0.015, 0.014, 0.041,
                              -0.166, 0.029, 0.278, 0.077, 0.208)

beta[c(1:10,14:nrow(beta)),6] <- c(0.031, -0.058, 0.230, 0.196, 0.342,
                              0.146, -0.035, -0.003, 0.096, 0.060,
                              0.013, 0.338, 0.201, 0.059, 0.135)

beta[c(1:10,14:nrow(beta)),7] <- c(0.129, -0.074, -0.246, -0.169, 0.070,
                              -0.375, -0.046, -0.385, -0.043, -0.015,
                              -0.052, 0.157, -0.024, -0.333, 0.621)

beta[c(1:10,14:nrow(beta)),8] <- c(-0.181, -0.061, 0.418, 0.413, 0.560,
                              0.366, 0.021, 0.109, -0.221, 0.045,
                              -0.283, 0.135, 0.038, 0.770, -0.369)

# Covariate generating distribution
add_cov <- list()

gen_L0 <- function(N) rbinom(N, 1, 0.4)                                                                       # koen
gen_A0 <- function(N, L0) pmin(rexp(N, 0.3) + 70, 100)                                                        # alder
add_cov[[1]] <- function(N) rbinom(N, 1, 0.16)                                                                # civst 1
add_cov[[2]] <- function(N) rbinom(N, 1, 0.11)                                                                # civst 2
add_cov[[3]] <- function(N) rbinom(N, 1, 0.6)                                                                 # civst 3
add_cov[[4]] <- function(N) rbinom(N, 1, 0.03)                                                                # civst 4
#add_cov[[6]] <- function(N) rbinom(N, 1, 0.01)                                                               # civst 5
add_cov[[5]] <- function(N) rpois(N, 0.25)                                                                    # n_diag_base
add_cov[[6]] <- function(N) pmax(rnorm(N, 2, 1), 0.2)                                                         # base_LDL
add_cov[[7]] <- function(N) rbinom(N, 1, 0.14)                                                                # A0
add_cov[[8]] <- function(N) rpois(N, 5)                                                                       # base_drugs

# Estimerede parametre, som vi har fået ved at fitte på event of interest og terminale events
nu <- c(2.4402339, 1.2381952, 0.8753749, 1.372, 0.757, 0.864, 0.623, 0.677)
eta <- c(4.352e-08, 6.689-04, 4.801-04, 5.525e-03, 1.697-02, 4.624e-02, 6.739e-02, 4.768e-02)


# Estimerede parametre, som vi har fået ved at fitte lm fits
eta <- c(10^-9, 0.001, 10^-9, 0.007, 0.032, 0.045, 0.063, 0.037)
nu <- c(1.853, 0.855, 0.763, 1.264, 0.743, 0.893, 0.711, 0.689)

# Simulating from simStatinData
data <- simStatinData(beta = beta, 
                      N = 2000, 
                      add_cov = add_cov, 
                      followup = 60, 
                      gen_A0 = gen_A0, 
                      gen_L0 = gen_L0,
                      eta = eta, 
                      nu = nu,
                      max_iter = 300,
                      lower = 10^(-30),
                      upper = 10^6,
                      cens = 0)

plotEventData(data[1:2000,])
data <- IntFormatData(data, N_cols = (n_cov + 4):(n_cov+n_proc+3))

#-------------------------------------------------------------------------------
# Fitting Models
#-------------------------------------------------------------------------------

# Fit models
# Models where the indidividuals are always at risk
vars <- c("L0", "A0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "OS", "L", "A", "LDL1", "LDL2")
form <- as.formula(paste("Surv(tstart, tstop, Delta == 0) ~", paste(vars, collapse = " + ")))
survfit1 <- coxph(form, data = data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 1) ~", paste(vars, collapse = " + ")))
survfit2 <- coxph(form, data = data)

form <- as.formula(paste("Surv(tstart, tstop, Delta == 2) ~", paste(vars, collapse = " + ")))
survfit3 <- coxph(form, data = data)

# Models for which the individual is not always at risk
form <- as.formula(paste("Surv(tstart, tstop, Delta == 3) ~", paste(vars, collapse = " + ")))
survfit4 <- coxph(form, data = data[OS <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 4) ~", paste(vars, collapse = " + ")))
survfit5 <- coxph(form, data = data[L <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 5) ~", paste(vars, collapse = " + ")))
survfit6 <- coxph(form, data = data[A <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 6) ~", paste(vars, collapse = " + ")))
survfit7 <- coxph(form, data = data[LDL1 <= 10])

form <- as.formula(paste("Surv(tstart, tstop, Delta == 7) ~", paste(vars, collapse = " + ")))
survfit8 <- coxph(form, data = data[LDL2 <= 10])

# Check estimations
beta_comp <- beta[c(1:10,14:nrow(beta)), ]

df <- rbind(cbind(beta_comp[,1],  confint(survfit1), 1), cbind(beta_comp[,2],  confint(survfit2), 2))
df <- rbind(df, cbind(beta_comp[,3],  confint(survfit3), 3))
df <- rbind(df, cbind(beta_comp[,4],  confint(survfit4), 4))
df <- rbind(df, cbind(beta_comp[,5],  confint(survfit5), 5))
df <- rbind(df, cbind(beta_comp[,6],  confint(survfit6), 6))
df <- rbind(df, cbind(beta_comp[,7],  confint(survfit5), 7))
df <- rbind(df, cbind(beta_comp[,8],  confint(survfit6), 8))


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
for(var in 4:(3+n_cov)){
  list_old_vars[(var - 3)] <- data[tstart == 0,..var]
}

names(list_old_vars) <- colnames(data[, 4:(3+n_cov)])

cox_fits <- list(
  "C" = survfit1,
  "D" = survfit2,
  "CVD" = survfit3,
  "LDL1" = survfit4,
  "LDL2" = survfit5,
  "OS" = survfit6,
  "A" = survfit7,
  "L" = survfit8
)


sim_data0 <- simEventCox(
  10^4,
  cox_fits,
  list_old_vars = list_old_vars,
  n_event_max = c(1, 1, 1, 10, 20, 10, 10, 10),
  term_events = c(1, 2, 3),
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
df <- rbind(cbind(beta[c(1:n_cov,(n_cov+4):(n_cov +n_proc)),1],  confint(survfit1), 1),
            cbind(beta[c(1:n_cov,(n_cov+4):(n_cov +n_proc)),2],  confint(survfit2), 2))

df <- rbind(df, cbind(beta[c(1:n_cov,(n_cov+4):(n_cov +n_proc)),3],  confint(survfit3), 3))
df <- rbind(df, cbind(beta[c(1:n_cov,(n_cov+5):(n_cov +n_proc)),4],  confint(survfit4), 4))
df <- rbind(df, cbind(beta[c(1:n_cov,(n_cov+4):(n_cov +n_proc)),5],  confint(survfit5), 5))
df <- rbind(df, cbind(beta[c(1:n_cov,(n_cov+4):(n_cov +n_proc)),6],  confint(survfit6), 6))


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

# An intervention multiplying the off statins process by alpha
alpha <- 2
intervention <- function(j, basehaz) if(j == 4) alpha * basehaz else basehaz

sim_data_int <- simEventCox(
  10^4,
  cox_fits,
  list_old_vars = list_old_vars,
  n_event_max = c(1, 1, 1, 1, 3, 3),
  term_events = c(1, 2, 3),
  intervention2 = intervention
)


ggarrange(plotEventData(sim_data0[ID %in% 1:200], title = "Non intervened") + xlim(c(0,1.5)), 
          plotEventData(sim_data_int[ID %in% 1:200], title = "Intervened")+ xlim(c(0,1.5)), ncol = 2)

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

alphas <- seq(0.5,2, by = 0.1)

#beta <-matrix(0, ncol = n_proc, nrow = n_cov+n_proc)
#beta[(n_cov+4), c(2,3)] <- 1

risk_alpha <- matrix(nrow = length(alphas), ncol = 2)
  
for(i in seq_along(alphas)){
  print(i)
  res_sim <- alphaSim(N = 1e5,
                      eta = rep(0.1,6),
                      nu = rep(1.1,6),
                      alpha = alphas[i],
                      tau = 5,
                      setting = "Statin",
                      cens = 0,
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


