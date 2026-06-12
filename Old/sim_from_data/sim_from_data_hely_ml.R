### sim_from_data_hely.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: May  9 2025 (20:17) 
## Version: 
## Last-Updated: May  9 2025 (20:17) 
##           By: Helene
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

set.seed(372)
beta <-  matrix(rnorm(3*2), nrow = 2)
data <- simCRdata(20000, beta =beta)
predict_data <- data.frame("L0" = 1, "A0" = 0.5)
true_surv_death_func <- function(t) exp( - 0.1  * t^(1.1) * exp(sum(beta[,2] * predict_data))) 
true_surv_t2d_func <- function(t) exp( - 0.1  * t^(1.1) * exp(sum(beta[,3] * predict_data)))

# Fit Cox models
cox_death <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = data)
cox_t2d <- coxph(Surv(Time, Delta == 2) ~ L0 + A0, data = data)
surv_death <- survfit(cox_death, newdata = predict_data, method = "breslow")
surv_t2d <- survfit(cox_t2d, newdata = predict_data, method = "breslow")

cox_fits <- list(death = cox_death, t2d = cox_t2d)
#new_data <- simEventDataCox(5000, cox_fits, L0_old = data$L0, term_events = c(1,2))


basehaz_list <- lapply(cox_fits, function(model) {
    bh_tmp <- basehaz(model, centered=FALSE)
    bh_tmp$dhazard <- diff(c(0,bh_tmp$hazard))
    return(bh_tmp)
})

N <- 2000

sim_data <- data.frame(
    L0 = sample(data$L0, N, replace = TRUE),  # sample with replacement from old data
    A0 = rbinom(N, 1, 0.5)
)

lp_alive <- lapply(cox_fits, function(model) predict(model, newdata = sim_data, type = "lp"))

# Jeg har kigget lidt på argumentet reference = "zero" i predict funktionen, men det gør vist ikke den store forskel

U <- matrix(-log(runif(N * num_events)), ncol = num_events) #-log(runif(N))
T_k <- rep(1, N)

num_events <- 2

event_times <- matrix(ncol = num_events, nrow = N)      # A matrix for event times

for (j in 1:num_events){
    event_times[,j] <- sapply(1:N, function(i) {
        lp_i_j <- lp_alive[[j]][i]
        basehaz_i_j <- basehaz_list[[j]][basehaz_list[[j]][, "time"]>T_k[i],]
        basehaz_i_j$haz <- cumsum(basehaz_i_j$dhazard*exp(lp_i_j))
        times_i_j <- basehaz_i_j$time[basehaz_i_j$haz<=U[i,j]]
        if(length(times_i_j)==0) {
            return(Inf) 
        } else {
            return(max(times_i_j))
        }
    })
}

sim_data <- setDT(sim_data)
sim_data[, Time:=#event_times[,2]
#event_times[,1]
apply(event_times, 1, min)
][, Delta:=1][Time==event_times[,2], Delta:=2]

new_data <- copy(sim_data)

new_data[, table(Delta)]
data[, table(Delta)]

# Fit Cox models
(cox_death2 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = new_data))
cox_death
(cox_t2d2<- coxph(Surv(Time, Delta == 2) ~ L0 + A0, data = new_data))
cox_t2d

surv_death2 <- survfit(cox_death2, newdata = predict_data, method = "breslow")
surv_t2d2 <- survfit(cox_t2d2, newdata = predict_data, method = "breslow")


ggplot()+
    geom_line(aes(x = surv_death$time, y = surv_death$surv), color = "red") +
    geom_line(aes(x = surv_death2$time, y = surv_death2$surv), color = "blue") +
    geom_function(fun = true_surv_death_func)

ggplot()+
    geom_line(aes(x = surv_t2d$time, y = surv_t2d$surv), color = "red") +
    geom_line(aes(x = surv_t2d2$time, y = surv_t2d2$surv), color = "blue") +
    geom_function(fun = true_surv_t2d_func)


######################################################################
### sim_from_data_hely.R ends here
