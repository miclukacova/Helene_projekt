#-------------------------------------------------------------------------------
# Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(simevent)

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
data <- simStatinData(beta = beta, N = 1000, add_cov = add_cov)
data <- IntFormatData(data, N_cols = 25:30)

#-------------------------------------------------------------------------------
# Fitting Models
#-------------------------------------------------------------------------------

# Fit models
data1 <- subset(data, select = -c(ID, Time, k, N0, N1, N2, N3))
survfit1 <- coxph(Surv(tstart, tstop, Delta == 0) ~ ., data = data1)
survfit2 <- coxph(Surv(tstart, tstop, Delta == 1) ~ ., data = data1)
survfit3 <- coxph(Surv(tstart, tstop, Delta == 2) ~ ., data = data1)
data4 <- subset(data, select = -c(ID, Time, k, N0, N1, N2, N3))
survfit4 <- coxph(Surv(tstart, tstop, Delta == 3) ~ . - ID -Time -k -N0 -N1 -N2 -N3, data = data[N3 == 0])
data5 <- subset(data, select = -c(ID, Time, k, N0, N1, N2, N4))
survfit5 <- coxph(Surv(tstart, tstop, Delta == 2) ~ . - ID -Time -k -N0 -N1 -N2 -N4, data = data[N4 <= 0])
data6 <- subset(data, select = -c(ID, Time, k, N0, N1, N2, N5))
survfit6 <- coxph(Surv(tstart, tstop, Delta == 3) ~ . - ID -Time -k -N0 -N1 -N2 -N5, data = data[N3 <= 0])
data[which(data$tstart >= data$tstop),]

#-------------------------------------------------------------------------------
# Simulating new data, respectively under intervention and without intervention
#-------------------------------------------------------------------------------


