########### SImulation using fitted RF ################

#-------------------------------------------------------------------------------
# Survival Data
#-------------------------------------------------------------------------------

set.seed(373)
beta = matrix(c(0, 0,0.5,1), ncol = 2, nrow = 2)
data <- simSurvData(N = 2*10^4, beta = beta, cens = 1)

# Fit model
RF_fit <- rfsrc(Surv(Time, Delta) ~ L0 + A0, data = data, ntree = 10, nodesize = 125)

new_data <- simEventRF(2*10^4, RF_fit, L0_old = data$L0, A0_old = data$A0,
                       n_event_max = 1, term_events = 1)


# Fit Cox models
cox1 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = new_data)
basehazz1 <- basehaz(cox1, centered = TRUE)

covars <- data.frame(L0 = c(0.5,0.5), A0 = c(1,0))

cox_term1 <- exp(predict(cox1, newdata=covars, type="lp"))
cum_int1 <- outer(basehazz1[['hazard']], cox_term1, "*")[1:9000,]

# Evaluate the true survival functions at these time points
true_chf_func1 <- function(t) 0.1  * t^(1.1) * exp(sum(covars[1,] * c(0.5,1)))
true_chf_func2 <- function(t) 0.1  * t^(1.1) * exp(sum(covars[2,] * c(0.5,1)))
time_points <- basehazz1[['time']][1:9000]
true_chf1 <- true_chf_func1(time_points)
true_chf2 <- true_chf_func2(time_points)


# Plot
par(mfrow = c(1,1), cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,6,1,1), mgp = c(4, 1, 0))
plot(time_points, cum_int1[,1], type="l", xlab="Time (Year)",
     ylab="Survival", col=1, lty=1, lwd=2)
lines(time_points, cum_int1[,2], col=2, lty=1, lwd=2)
# Add lines for the true survival functions
lines(time_points, true_chf1, col=1, lty=3, lwd=2)
lines(time_points, true_chf2, col=2, lty=3, lwd=2)
legend("topright",
       legend=c("A0=1 (Predicted)", "A0=0 (Predicted)", "A0=1 (True)", "A0=0 (True)"),
       col=c(1, 2, 1, 2),
       lty=c(1, 2, 3, 4),
       cex=0.75,
       lwd=2)

# Måske er der noget som går galt her?

#-------------------------------------------------------------------------------
#### Competing Risk
#-------------------------------------------------------------------------------

# Simulate data
set.seed(373)
beta = matrix(c(0.5,-1,-0.5,0.5,0,0.5), ncol = 3, nrow = 2)
data <- simCRdata(N = 10^4, beta = beta)

# Fit RSF model
#tune.nodesize(Surv(Time,Delta) ~ L0 + A0, data = data)
RF_fit <- rfsrc(Surv(Time, Delta) ~ L0 + A0, data = data, ntree = 10, nodesize = 150)

# New Data
new_data <- simEventRF(10^4, RF_fit, L0_old = data$L0, A0_old = data$A0,
           n_event_max = c(1,1), term_events = c(1,2))

# Fit Cox models
cox1 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = new_data)
cox2 <- coxph(Surv(Time, Delta == 2) ~ L0 + A0, data = new_data)
basehazz1 <- basehaz(cox1, centered = TRUE)
basehazz2 <- basehaz(cox2, centered = TRUE)
cox_term1 <- exp(predict(cox1, newdata=newdata, type="lp"))
cox_term2 <- exp(predict(cox2, newdata=newdata, type="lp"))
cum_int1 <- outer(basehazz1[['hazard']], cox_term1, "*")[1:9000,]
cum_int2 <- outer(basehazz2[['hazard']], cox_term2, "*")[1:9000,]

# True CHF
true_chf_func1 <- function(t, cov)  0.1  * t^(1.1) * exp(sum(cov * c(-0.5,0.5)))
true_chf_func2 <- function(t, cov)  0.1  * t^(1.1) * exp(sum(cov * c(0,0.5)))

# Evaluate the true survival functions at time points
times <- basehazz1[['time']][1:9000]
true_chf11 <- true_chf_func1(times, newdata[1,])
true_chf12 <- true_chf_func1(times, newdata[2,])
true_chf21 <- true_chf_func2(times, newdata[1,])
true_chf22 <- true_chf_func2(times, newdata[2,])


# Plot
par(mfrow = c(1,2), cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,6,1,1), mgp = c(4, 1, 0))
# Process 1
plot(times, true_chf11, type="l", xlab="Time",
     ylab="Process 1 CHF", col=1, lty=1, lwd=2)
lines(times, true_chf12, col=1, lty=2, lwd=2)
lines(times, cum_int1[,1], col=2, lty=1, lwd=2)
lines(times, cum_int1[,2], col=2, lty=2, lwd=2)

# Process 2
plot(times,true_chf21, type="l", xlab="Time",
     ylab="Process 2 CHF", col=1, lty=1, lwd=2)
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
