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
# There are 8 processes: Censoring, Death, CVD
n_cov <- 10
n_proc <- 3

# Regression coefficients

# Creating beta matrix
beta <- matrix(0, nrow = n_cov + n_proc, ncol = n_proc)

rownames(beta) <- c(
  "L0", "A0",
  "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8",
  "C", "D", "CVD")
colnames(beta) <- c("C", "D", "CVD")

#beta[1:10,1] <- c(0, 0.001, 0.021, 0.005, -0.008, -0.081, -0.008, 0.009, -0.023, 0.014)

beta[1:10,2] <- c(-0.492, -0.008, 0.451, 0.116, 0.657, 0.577, 0.049, -0.093, 2.650, 0.121)

beta[1:10,3] <- c(-0.614, -0.091, 0.731, 0.226, 0.747, 0, 0.305, -0.065, 0.197, 0.079)


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
eta <- c(1, 0.0000412, 0.0001225)
nu <- c(1, 1.197604, 1.014444)
at_risk = function(events, covariates) return(c(0,1,1))

# Estimerede parametre, som vi har fået ved at fitte lm fits
#eta <- rep(0.001,8)
#nu <- rep(1.001,8)

# Simulating from simStatinData
data <- simEventData(beta = beta, 
                     N = 2000, 
                     term_deltas = c(0, 1, 2),
                     add_cov = add_cov, 
                     max_cens = 60, 
                     gen_A0 = gen_A0, 
                     gen_L0 = gen_L0,
                     eta = eta, 
                     nu = nu,
                     max_iter = 300,
                     lower = 10^(-40),
                     upper = 10^8,
                     at_risk = at_risk)

plotEventData(data[1:2000,])
#data <- IntFormatData(data, N_cols = (n_cov + 4):(n_cov+n_proc+3))

#-------------------------------------------------------------------------------
# Fitting Models
#-------------------------------------------------------------------------------

# Fit models
# Models where the indidividuals are always at risk
vars <- c("L0", "A0", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")
form <- as.formula(paste("Surv(Time, Delta == 0) ~", paste(vars, collapse = " + ")))
survfit1 <- coxph(form, data = data)

form <- as.formula(paste("Surv(Time, Delta == 1) ~", paste(vars, collapse = " + ")))
survfit2 <- coxph(form, data = data)

form <- as.formula(paste("Surv(Time, Delta == 2) ~", paste(vars, collapse = " + ")))
survfit3 <- coxph(form, data = data)

# Check estimations

df <- rbind(cbind(beta[1:10,1],  confint(survfit1), 1), cbind(beta[1:10,2],  confint(survfit2), 2))
df <- rbind(df, cbind(beta[1:10,3],  confint(survfit3), 3))

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
  "CVD" = survfit3
)


sim_data0 <- simEventCox(
  10^4,
  cox_fits,
  list_old_vars = list_old_vars,
  n_event_max = c(1, 1, 1),
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

form <- as.formula(paste("Surv(Time, Delta == 0) ~", paste(vars, collapse = " + ")))
survfit1 <- coxph(form, data = sim_data)

form <- as.formula(paste("Surv(Time, Delta == 1) ~", paste(vars, collapse = " + ")))
survfit2 <- coxph(form, data = sim_data)

form <- as.formula(paste("Surv(Time, Delta == 2) ~", paste(vars, collapse = " + ")))
survfit3 <- coxph(form, data = sim_data)

# Check estimations
df <- rbind(cbind(beta[1:10,1],  confint(survfit1), 1),
            cbind(beta[1:10,2],  confint(survfit2), 2))
df <- rbind(df, cbind(beta[1:10,3],  confint(survfit3), 3))

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

