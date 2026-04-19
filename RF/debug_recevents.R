# Test of simulation functions

library(simevent)
library(ggplot2)
library(ggpubr)
library(data.table)

N <- 10^5

# It works fine on the whole data, and also on the splitted data

# In simDisease

# The observed data
data_obs <- simDisease(N)

data1 <- data_obs[order(Time), .SD[1], by = ID]
data2 <- data_obs[ID %in% data1[!(Delta %in% c(0,1)), ID], ]
data2 <- data2[order(Time), .SD[1:2], by = ID]
data2[, oldTime := min(Time), by = ID]
data2[, newTime := Time - oldTime, by = ID]
data2 <- data2[order(Time), .SD[2], by = ID]
setkey(data2, ID)

# Fit some Cox models
cox_cens <- survival::coxph(survival::Surv(Time, Delta == 0) ~ L0 + A0, data = data2)
cox_death <- survival::coxph(survival::Surv(Time, Delta == 1) ~ L0 + A0, data = data2)
#cox_Disease <- survival::coxph(survival::Surv(Time, Delta == 2) ~ L0 + A0, data = data_obs)

# Then simulate new data:
cox_fits <- list("C" = cox_cens, "D" = cox_death) #, "L" = cox_Disease)
list_old_vars <- list("L0" = data_obs$L0, "A0" = data_obs$A0)
new_data <- simEventCox(N, cox_fits = cox_fits, list_old_vars = list_old_vars,
                        term_events = c(1,2), n_event_max = c(1,1))

# Fit new Cox models
cox_cens2 <- survival::coxph(survival::Surv(Time, Delta == 0) ~ L0 + A0, data = new_data)
cox_death2 <- survival::coxph(survival::Surv(Time, Delta == 1) ~ L0 + A0, data = new_data)
#cox_Disease2 <- survival::coxph(survival::Surv(Time, Delta == 2) ~ L0 + A0, data = new_data)

# Compute confidence intervals
ci_cens <- confint(cox_cens2)
coef_cens <- cox_cens$coefficients

ci_death <- confint(cox_death2)
coef_death <- cox_death$coefficients

#ci_Disease <- confint(cox_Disease2)
#coef_Disease <- cox_Disease$coefficients

# Compare confidence intervals and true values
all(coef_cens >= ci_cens[, 1] & coef_cens <= ci_cens[, 2])
all(coef_death >= ci_death[, 1] & coef_death <= ci_death[, 2])
#all(coef_Disease >= ci_Disease[, 1] & coef_Disease <= ci_Disease[, 2])


ggarrange(plotEventData(new_data[1:2000,], title = "Sim Data"),
          plotEventData(data2[1:2000,], title = "Original Data"), ncol = 2)

new_data[, table(Delta)]/nrow(new_data)
data2[, table(Delta)]/nrow(data2)
new_data[, summary(Time)]
data2[, summary(Time)]


# In simDisease

# The observed data
data_obs <- simDisease(N)
data_obs <- IntFormatData(data_obs, N_cols = 6)

# Fit some Cox models
cox_cens <- survival::coxph(survival::Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + L, data = data_obs)
cox_death <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = data_obs)
cox_Disease <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = data_obs[L == 0])

# Then simulate new data:
cox_fits <- list("C" = cox_cens, "D" = cox_death, "L" = cox_Disease)
list_old_vars <- list("L0" = data_obs$L0, "A0" = data_obs$A0)
new_data <- simEventCox(N, cox_fits = cox_fits, list_old_vars = list_old_vars,
                        term_events = c(1,2), n_event_max = c(1,1,1))
new_data <- IntFormatData(new_data, N_cols = 6:8)

# Fit new Cox models
cox_cens2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + L, data = new_data)
cox_death2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = new_data)
cox_Disease2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = new_data[L == 0])

# Compute confidence intervals
ci_cens <- confint(cox_cens2)
coef_cens <- cox_cens$coefficients

ci_death <- confint(cox_death2)
coef_death <- cox_death$coefficients

ci_Disease <- confint(cox_Disease2)
coef_Disease <- cox_Disease$coefficients

# Compare confidence intervals and true values
all(coef_cens >= ci_cens[, 1] & coef_cens <= ci_cens[, 2])
all(coef_death >= ci_death[, 1] & coef_death <= ci_death[, 2])
all(coef_Disease >= ci_Disease[, 1] & coef_Disease <= ci_Disease[, 2])


ggarrange(plotEventData(new_data[1:2000,], title = "Sim Data"),
          plotEventData(data_obs[1:2000,], title = "Original Data"), ncol = 2)

new_data[, table(Delta)]/nrow(new_data)
data_obs[, table(Delta)]/nrow(data_obs)
new_data[, summary(Time)]
data_obs[, summary(Time)]


# In simTreatment

# The observed data
set.seed(373)
data_obs <- simTreatment(N, cens = 0)
data_obs <- IntFormatData(data_obs, N_cols = 5:6)

# Fit some Cox models
cox_D <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + L + A, data = data_obs)
cox_A <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + L, data = data_obs[A == 0])
cox_L <- survival::coxph(survival::Surv(tstart, tstop, Delta == 3) ~ L0 + A, data = data_obs[L == 0])

# Then simulate new data:
cox_fits <- list("D" = cox_D, "A" = cox_A, "L" = cox_L)
list_old_vars <- list("L0" = data_obs$L0)
new_data <- simEventCox(N, cox_fits = cox_fits, list_old_vars = list_old_vars,
                        term_events = c(1), n_event_max = c(1,1,1))
new_data <- IntFormatData(new_data, N_cols = 5:7)

# Fit new Cox models
cox_D2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 0) ~ L0 + L + A, data = new_data)
cox_A2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + L, data = new_data[A == 0])
cox_L2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A, data = new_data[L == 0])

# Compute confidence intervals
ci_A <- confint(cox_A2)
coef_A <- cox_A$coefficients

ci_L <- confint(cox_L2)
coef_L <- cox_L$coefficients

ci_D <- confint(cox_D2)
coef_D <- cox_D$coefficients

# Compare confidence intervals and true values
all(coef_L >= ci_L[, 1] & coef_L <= ci_L[, 2])
all(coef_A >= ci_A[, 1] & coef_A <= ci_A[, 2])
all(coef_D >= ci_D[, 1] & coef_D <= ci_D[, 2])

ggarrange(plotEventData(new_data[1:2000,], title = "Sim Data"),
          plotEventData(data_obs[1:2000,], title = "Original Data"), ncol = 2)

new_data[, table(Delta)]/nrow(new_data)
data_obs[, table(Delta)]/nrow(data_obs)
new_data[, summary(Time)]
data_obs[, summary(Time)]


# DET HER VIRKER RET GODT! PRØV DET SAMME MED SIMEVENTOBJ

# In simDisease

# The observed data
data_obs <- simDisease(N)
data_obs <- IntFormatData(data_obs, N_cols = 6)

# Fit some Cox models
cox_cens <- survival::coxph(survival::Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + L, data = data_obs)
cox_death <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = data_obs)
cox_Disease <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = data_obs[L == 0])

# Equip with predict2 method
predict2 <- function(obj, ...) {
  UseMethod("predict2")
}
predict2.simevent <- function(obj, sim_data){
  # Base hazards
  basehazz_list <- lapply(obj, function(model) basehaz(model, centered = FALSE))
  
  # Individual specific term
  cox_term <- lapply(obj, function(model)
    exp(stats::predict(model, newdata = sim_data, type="lp", reference = "zero")))
  
  # Chf
  chf_list <- lapply(1:length(obj), function(j) cox_term[[j]] %*% t(basehazz_list[[j]][["hazard"]]))
  chf <- array(dim = c(c(dim(chf_list[[1]])), length(obj)))
  for(j in 1:length(obj)) chf[,,j] <- chf_list[[j]]
  
  # Return list
  list(time = basehazz_list[[1]][["time"]], chf = chf)
}


# Then simulate new data:
cox_fits <- list("C" = cox_cens, "D" = cox_death, "L" = cox_Disease)
class(cox_fits) <- "simevent"
list_old_vars <- list("L0" = data_obs$L0, "A0" = data_obs$A0)
new_data <- simEventObj(N, obj = cox_fits, list_old_vars = list_old_vars, event_names = c("C", "D", "L"))


## Du er kommet her til....
  
new_data <- IntFormatData(new_data, N_cols = 6:8)

# Fit new Cox models
cox_cens2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + L, data = new_data)
cox_death2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = new_data)
cox_Disease2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = new_data[L == 0])

# Compute confidence intervals
ci_cens <- confint(cox_cens2)
coef_cens <- cox_cens$coefficients

ci_death <- confint(cox_death2)
coef_death <- cox_death$coefficients

ci_Disease <- confint(cox_Disease2)
coef_Disease <- cox_Disease$coefficients

# Compare confidence intervals and true values
all(coef_cens >= ci_cens[, 1] & coef_cens <= ci_cens[, 2])
all(coef_death >= ci_death[, 1] & coef_death <= ci_death[, 2])
all(coef_Disease >= ci_Disease[, 1] & coef_Disease <= ci_Disease[, 2])


ggarrange(plotEventData(new_data[1:2000,], title = "Sim Data"),
          plotEventData(data_obs[1:2000,], title = "Original Data"), ncol = 2)

new_data[, table(Delta)]/nrow(new_data)
data_obs[, table(Delta)]/nrow(data_obs)
new_data[, summary(Time)]
data_obs[, summary(Time)]

# In simTreatment

# The observed data
set.seed(373)
data_obs <- simTreatment(N, cens = 0)
data_obs <- IntFormatData(data_obs, N_cols = 5:6)

# Fit some Cox models
cox_D <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + L + A, data = data_obs)
cox_A <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + L, data = data_obs[A == 0])
cox_L <- survival::coxph(survival::Surv(tstart, tstop, Delta == 3) ~ L0 + A, data = data_obs[L == 0])

# Then simulate new data:
cox_fits <- list("D" = cox_D, "A" = cox_A, "L" = cox_L)
class(cox_fits) <- "simevent"
list_old_vars <- list("L0" = data_obs$L0)
new_data <- simEventObj(N, obj = cox_fits, list_old_vars = list_old_vars, event_names("D", "A", "L"))
new_data <- IntFormatData(new_data, N_cols = 5:7)

# Fit new Cox models
cox_D2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 0) ~ L0 + L + A, data = new_data)
cox_A2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + L, data = new_data[A == 0])
cox_L2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A, data = new_data[L == 0])

# Compute confidence intervals
ci_A <- confint(cox_A2)
coef_A <- cox_A$coefficients

ci_L <- confint(cox_L2)
coef_L <- cox_L$coefficients

ci_D <- confint(cox_D2)
coef_D <- cox_D$coefficients

# Compare confidence intervals and true values
all(coef_L >= ci_L[, 1] & coef_L <= ci_L[, 2])
all(coef_A >= ci_A[, 1] & coef_A <= ci_A[, 2])
all(coef_D >= ci_D[, 1] & coef_D <= ci_D[, 2])

ggarrange(plotEventData(new_data[1:2000,], title = "Sim Data"),
          plotEventData(data_obs[1:2000,], title = "Original Data"), ncol = 2)

new_data[, table(Delta)]/nrow(new_data)
data_obs[, table(Delta)]/nrow(data_obs)
new_data[, summary(Time)]
data_obs[, summary(Time)]




