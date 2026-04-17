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
# There are 6 processes: Censoring, Death, CVD, Off Statins, Treatment, Disease

n_cov <- 7
n_proc <- 6

# Regression coefficients

beta <- matrix(0, nrow = n_cov + n_proc, ncol = n_proc)
beta[,1] <- c(0.07, -0.03, 0, )
beta[,2] <- c(0, 0, 0.5, 0, 0.1, -0.1, 0.1, 0,0,0, -0.15, 1.1, 0.6)
beta[,3] <- c(0, 0.35, 0.5, 0, 0.05, 0, 0, 0,0,0, -0.05, 0.5, 0)

# Covariate generating distribution
add_cov <- list()
gen_A0 <- function(N, L0) rbinom(N, 1, 0.14)                                                                  # A0
gen_L0 <- function(N) rpois(N, 0.1)                                                                           # base_diags
add_cov[[1]] <- function(N) rbinom(N, 1, 0.4)                                                                 # Sex
add_cov[[2]] <- function(N) sample(x = 1:4, size = N, prob = c(0.6, 0.1, 0.1, 0.2), replace = TRUE)           # civst
add_cov[[3]] <- function(N) rexp(N, 0.3) + 70                                                                 # Age
add_cov[[4]] <- function(N) pmax(rnorm(N, 2, 1), 0.5)                                                         # base_LDL
add_cov[[5]] <- function(N) rpois(N, 5)                                                                       # base_drugs

# Simulating from simStatinData
data <- simStatinData(beta = beta, N = 5*10^4, add_cov = add_cov, followup = Inf, 
                      gen_A0 = gen_A0, eta = rep(0.01, 6), nu = rep(1.02,6))
data <- IntFormatData(data, N_cols = (n_cov + 4):(n_cov+n_proc+3))

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

vars <- setdiff(names(data), c("ID", "Time", "k", "C", "D", "CVD", "Delta", "tstart", "tstop"))
form <- as.formula(paste("Surv(tstart, tstop, Delta == 5) ~", paste(vars, collapse = " + ")))

survfit6 <- coxph(form, data = data[L <= 3])

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
# Simulating new data without intervention
#-------------------------------------------------------------------------------

list_old_vars <- list()
for(var in 4:(3+n_cov)){
  # Hvad er den rigtige måde at trække på her?
  #list_old_vars[(var - 3)] <- data[,..var]
  list_old_vars[(var - 3)] <- data[tstart == 0,..var]
}

names(list_old_vars) <- colnames(data[, 4:(3+n_cov)])

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


