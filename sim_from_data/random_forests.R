#### Cooking with the random forest fits ####

library(randomForestSRC)
library(devtools)
load_all()

#-------------------------------------------------------------------------------
# Simulate Survival Data
#-------------------------------------------------------------------------------

set.seed(373)
beta = matrix(c(0, 0,0.5,-0.5), ncol = 2, nrow = 2)
data <- simSurvData(N = 1000, beta = beta)

# Fit model
tune.nodesize(Surv(Time,Delta) ~ L0 + A0, data = data)
obj <- rfsrc(Surv(Time, Delta == 1) ~ L0 + A0, data = data, ntree = 10, nodesize = 125)            
#obj <- rfsrc.fast(Surv(Time, Delta == 1) ~ L0 + A0, data = data)   

# plot tree number 3
plot(get.tree(obj, 3))

## plot survival curves for first 10 individuals -- direct way
matplot(obj$time.interest, 100 * t(obj$survival.oob[1:10, ]),
        xlab = "Time", ylab = "Survival", type = "l", lty = 1)

newdata <- data.frame(L0 = c(0.5,0.5), A0 = c(1,0))
y.pred <- predict(obj, newdata = newdata)
true_surv_death_func1 <- function(t) exp( - 0.1  * t^(1.1) * exp(sum(newdata[1,] * c(0.5,-0.5))))
true_surv_death_func2 <- function(t) exp( - 0.1  * t^(1.1) * exp(sum(newdata[2,] * c(0.5,-0.5))))

# Evaluate the true survival functions at these time points
time_points <- round(y.pred$time.interest, 2)
true_surv1 <- true_surv_death_func1(time_points)
true_surv2 <- true_surv_death_func2(time_points)


# Plot
par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,6,1,1), mgp = c(4, 1, 0))
plot(round(y.pred$time.interest,2),y.pred$survival[1,], type="l", xlab="Time (Year)",
     ylab="Survival", col=1, lty=1, lwd=2)
lines(round(y.pred$time.interest,2), y.pred$survival[2,], col=2, lty=2, lwd=2)
# Add lines for the true survival functions
lines(time_points, true_surv1, col=1, lty=3, lwd=2)  # Dashed line for func1
lines(time_points, true_surv2, col=2, lty=4, lwd=2)  # Dot-dashed line for func2
legend("topright", 
       legend=c("A0=1 (Predicted)", "A0=0 (Predicted)", "A0=1 (True)", "A0=0 (True)"), 
       col=c(1, 2, 1, 2), 
       lty=c(1, 2, 3, 4), 
       cex=2, 
       lwd=2)

#-------------------------------------------------------------------------------
#### Competing Risk
#-------------------------------------------------------------------------------

# Simulate data
set.seed(373)
beta = matrix(c(0.5,-1,-0.5,0.5,0,0.5), ncol = 3, nrow = 2)
data <- simCRdata(N = 10^4, beta = beta)

# Fit RSF model
#tune.nodesize(Surv(Time,Delta) ~ L0 + A0, data = data)
obj0 <- rfsrc(Surv(Time, Delta) ~ L0 + A0, data = data, ntree = 10, nodesize = 150) 

# New Data
newdata <- data.frame(L0 = c(0.5,0.5), A0 = c(1,0))
y.pred <- predict(obj0, newdata = newdata)

# True CHF
true_chf_func1 <- function(t, cov)  0.1  * t^(1.1) * exp(sum(cov * c(-0.5,0.5)))
true_chf_func2 <- function(t, cov)  0.1  * t^(1.1) * exp(sum(cov * c(0,0.5)))

# Evaluate the true survival functions at time points
time_points <- round(y.pred$time.interest, 2)
true_chf11 <- true_chf_func1(time_points, newdata[1,])
true_chf12 <- true_chf_func1(time_points, newdata[2,])
true_chf21 <- true_chf_func2(time_points, newdata[1,])
true_chf22 <- true_chf_func2(time_points, newdata[2,])

# Fit Cox models
cox1 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = data)
cox2 <- coxph(Surv(Time, Delta == 2) ~ L0 + A0, data = data)
basehazz1 <- basehaz(cox1, centered = TRUE)
basehazz2 <- basehaz(cox2, centered = TRUE)
cox_term1 <- exp(predict(cox1, newdata=newdata, type="lp"))
cox_term2 <- exp(predict(cox2, newdata=newdata, type="lp"))
cum_int1 <- outer(basehazz1[['hazard']], cox_term1, "*")
cum_int2 <- outer(basehazz2[['hazard']], cox_term2, "*")


# Plot
par(mfrow = c(1,2), cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,6,1,1), mgp = c(4, 1, 0))
# Process 1
plot(y.pred$time.interest,y.pred$chf[1,,1], type="l", xlab="Time",
     ylab="Process 1 CHF", col=1, lty=1, lwd=2)
lines(y.pred$time.interest, y.pred$chf[2,,1], col=1, lty=2, lwd=2)
# Add lines for the true survival functions
lines(time_points, true_chf11, col=2, lty=1, lwd=2)  
lines(time_points, true_chf12, col=2, lty=2, lwd=2)  
lines(basehazz1$time, cum_int1[,1], col=3, lty=1, lwd=2)
lines(basehazz2$time, cum_int1[,2], col=3, lty=2, lwd=2)

# Process 2
plot(y.pred$time.interest,y.pred$chf[1,,2], type="l", xlab="Time",
     ylab="Process 2 CHF", col=1, lty=1, lwd=2)
lines(y.pred$time.interest, y.pred$chf[2,,2], col=1, lty=2, lwd=2)
# Add lines for the true survival functions
lines(time_points, true_chf21, col=2, lty=1, lwd=2)  
lines(time_points, true_chf22, col=2, lty=2, lwd=2) 
lines(basehazz2$time, cum_int2[,1], col=3, lty=1, lwd=2)
lines(basehazz2$time, cum_int2[,2], col=3, lty=2, lwd=2)

legend("topright", 
       legend = c("RF CHF A0 = 1", "RF CHF A0 = 0", "True CHF A0 = 1", "True CHF A0 = 0",
                  "Cox CHF A0 = 1", "Cox CHF A0 = 0"),
       col = c(1, 1, 2, 2, 3, 3),
       lty = c(1, 2, 1, 2, 1, 2),
       lwd = 2,
       cex = 0.75)


# Den prædikterer dårligere:(((


#-------------------------------------------------------------------------------
# General Event Data
#-------------------------------------------------------------------------------

# recforest?

set.seed(373)
data <- simDisease(N = 5000, beta_L0_L = 1, beta_L0_D = 1, beta_L_D = 0.5, beta_A0_L = 1)
data_int <- IntFormatData(data, N_cols = 6)

# Fit models
rf1 <- rfsrc(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = data_int, ntree = 10, nodesize = 150)
rf2 <- rfsrc(Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = data_int, ntree = 10, nodesize = 150)

plot(rf1)
plot(rf2)



