# ------------------------------------------------------------------------------
# Object based on RF
beta = matrix(c(0.5,-1,-0.5,0.5,0,0.5), ncol = 3, nrow = 2)
data <- simCRdata(N = 5000, beta = beta)
list_old_vars <- list(L0 = data$L0, A0 = data$A0)

RF_fit <- randomForestSRC::rfsrc(Surv(Time, Delta) ~ L0 + A0, data = data)

predict2 <- function(obj, ...) {
  UseMethod("predict2")
}
predict2.rfsrc <- function(obj, sim_data, ...){
  preds <- stats::predict(RF_fit, sim_data)
  if(is.na(dim(preds$chf)[3])){
    chf <- array(dim = c(dim(preds$chf), 1))
    chf[,,1] <- preds$chf
  } else {
    chf <- preds$chf
  }
  list(time = preds$time.interest, chf = chf)
}

# New data
new_data <- simEventObj(5000, RF_fit, list_old_vars = list_old_vars)

ggarrange(plotEventData(new_data, title = "Sim Data"),
          plotEventData(data, title = "Original Data"), ncol = 2)

new_data[, table(Delta)]/nrow(new_data)
data_obs[, table(Delta)]/nrow(data_obs)
new_data[, summary(Time)]
data_obs[, summary(Time)]


# ------------------------------------------------------------------------------
# Object based on Cox
# ------------------------------------------------------------------------------
# Simulate CR data
# ------------------------------------------------------------------------------
set.seed(373)
beta = matrix(c(0.5,-1,-0.5,0.5,0,0.5), ncol = 3, nrow = 2)
data <- simCRdata(N = 3000, beta = beta)
list_old_vars <- list(L0 = data$L0, A0 = data$A0)

# Fit Cox Models
cox1 <- survival::coxph(survival::Surv(Time, Delta == 1) ~ L0 + A0, data = data)
cox2 <- survival::coxph(survival::Surv(Time, Delta == 2)~ L0 + A0, data = data)

# Create Object
cox_fits <- list("D" = cox1, "L" = cox2)
class(cox_fits) <- "simevent"

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

# Simulate new data
new_data <- simEventObj(10^4, cox_fits, list_old_vars = list_old_vars)

# Check whether new data corresponds to the old data
covars <- data.frame(L0 = c(0.5,0.5), A0 = c(1,0))
coxfit1 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = new_data)
coxfit2 <- coxph(Surv(Time, Delta == 2) ~ L0 + A0, data = new_data)
basehazz1 <- basehaz(coxfit1, centered = TRUE)
basehazz2 <- basehaz(coxfit2, centered = TRUE)
cox_term1 <- exp(predict(coxfit1, newdata=covars, type="lp"))
cox_term2 <- exp(predict(coxfit2, newdata=covars, type="lp"))
cum_int1 <- outer(basehazz1[['hazard']], cox_term1, "*")
cum_int2 <- outer(basehazz2[['hazard']], cox_term2, "*")
true_chf_func1 <- function(t, cov)  0.1  * t^(1.1) * exp(sum(cov * c(-0.5,0.5)))
true_chf_func2 <- function(t, cov)  0.1  * t^(1.1) * exp(sum(cov * c(0,0.5)))

# Evaluate the true survival functions at time points
times <- basehazz1[['time']]
true_chf11 <- true_chf_func1(times, covars[1,])
true_chf12 <- true_chf_func1(times, covars[2,])
true_chf21 <- true_chf_func2(times, covars[1,])
true_chf22 <- true_chf_func2(times, covars[2,])

par(mfrow = c(1,2), cex.axis = 1.0, cex.lab = 1.0, cex.main = 1.0, mar = c(6.0,6,1,1), mgp = c(4, 1, 0))
# Process 1
plot(times, true_chf11, type="l", xlab="Time",
     ylab="Process 1 CHF", col=1, lty=1, lwd=2,
     xlim = c(0,quantile(new_data$Time, 0.95)),
     ylim = c(0,2))
lines(times, true_chf12, col=1, lty=2, lwd=2)
lines(times, cum_int1[,1], col=2, lty=1, lwd=2)
lines(times, cum_int1[,2], col=2, lty=2, lwd=2)

# Process 2
plot(times,true_chf21, type="l", xlab="Time",
     ylab="Process 2 CHF", col=1, lty=1, lwd=2,
     xlim = c(0,quantile(new_data$Time, 0.95)),
     ylim = c(0,3))
lines(times, true_chf22, col=1, lty=2, lwd=2)
lines(times, cum_int2[,1], col=2, lty=1, lwd=2)
lines(times, cum_int2[,2], col=2, lty=2, lwd=2)

legend("topright",
       legend = c("True CHF A0 = 1", "True CHF A0 = 0",
                  "Cox CHF A0 = 1", "Cox CHF A0 = 0"),
       col = c(1, 1, 2, 2),
       lty = c(1, 2, 1, 2),
       lwd = 2,
       cex = 0.75)

# ------------------------------------------------------------------------------
# Survival Data
# ------------------------------------------------------------------------------

beta = matrix(c(1, 0,0.5,1), ncol = 2, nrow = 2)
data <- simSurvData(N = 3000, beta = beta, cens = 0)
list_old_vars <- list(L0 = data$L0, A0 = data$A0)

# Fit Cox Model
cox1 <- survival::coxph(survival::Surv(Time, Delta == 1) ~ L0 + A0, data = data)

# Create Object
cox_fits <- list("D" = cox1)
class(cox_fits) <- "simevent"

new_data <- simEventObj(2*10^4, cox_fits, list_old_vars = list_old_vars)
covars <- data.frame(L0 = c(0.5,0.5), A0 = c(1,0))

cox1 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = new_data)
basehazz1 <- basehaz(cox1, centered = TRUE)
cox_term1 <- exp(predict(cox1, newdata=covars, type="lp"))
cum_int1 <- outer(basehazz1[['hazard']], cox_term1, "*")

true_chf_func1 <- function(t) 0.1  * t^(1.1) * exp(sum(covars[1,] * c(0.5,1)))
true_chf_func2 <- function(t) 0.1  * t^(1.1) * exp(sum(covars[2,] * c(0.5,1)))
time_points <- basehazz1[['time']]
true_chf1 <- true_chf_func1(time_points)
true_chf2 <- true_chf_func2(time_points)

par(mfrow = c(1,1), cex.axis = 1.0, cex.lab = 1.0, cex.main = 1.0, mar = c(6.0,6,1,1), mgp = c(4, 1, 0))
plot(time_points, cum_int1[,1], type="l", xlab="Time (Year)",
     ylab="CHF", col=1, lty=1, lwd=2, xlim = c(0,quantile(new_data$Time, 0.9)),
     ylim = c(0,5))
lines(time_points, cum_int1[,2], col=2, lty=1, lwd=2)
# Add lines for the true survival functions
lines(time_points, true_chf1, col=1, lty=3, lwd=2)
lines(time_points, true_chf2, col=2, lty=3, lwd=2)
legend("topright",
       legend=c("A0=1 (Predicted)", "A0=0 (Predicted)", "A0=1 (True)", "A0=0 (True)"),
       col=c(1, 2, 1, 2),
       lty=c(1, 1, 3, 3),
       cex=0.75,
       lwd=2)














