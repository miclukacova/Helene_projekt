# Summary Function for proportions

summary_sim_event <- function(data) {
  N <- max(data$ID)
  prop_op <- nrow(data[Delta == 0])/N
  prop_death <- nrow(data[Delta == 1])/N
  prop_cens <- nrow(data[Delta == 2])/N
  prop_fu_cens <- sum(is.na(data$Delta))/N
  prop_cov <-  nrow(data[Delta == 3])/N

  mean_T <- mean(data$Time)
  max_T <- max(data$Time)

  return(c("Prop Op" = prop_op, "Prop Death" = prop_death, "Prop Cens" = prop_cens,
           "Prop FU Cens" = prop_fu_cens, "Prop Cov" = prop_cov, "Mean Time" = mean_T,
           "Max Time" = max_T))
}

# Function for fitting models, with and without misspecifications

fit_mod_set1 <- function(data, timefix = TRUE){
  trans_data <- trans_int_data(data)

  trans_data[, at_risk_oper := as.numeric(A == 0)]
  trans_data[, at_risk_cov := as.numeric(L == 0)]

  results <- matrix(nrow = 10, ncol = 3)
  colnames(results) <- c("Fully Specified", "-A", "-L")
  rownames(results) <- c("L0", "L", "L0", "A", "L", "L0", "A", "L", "L0", "A")

  # True models:
  survfit_oper <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + L,
                        data = trans_data[at_risk_oper == 1], cluster = ID,
                        timefix = timefix)
  survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A + L,
                         data = trans_data, cluster = ID,
                         timefix = timefix)
  survfit_cens <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A + L,
                        data = trans_data, cluster = ID,
                        timefix = timefix)
  survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + A,
                       data = trans_data[at_risk_cov == 1], cluster = ID,
                       timefix = timefix)

  results[,1] <- c(survfit_oper$coefficients, survfit_death$coefficients,
                   survfit_cens$coefficients, survfit_cov$coefficients)

  # Cox models only with time varying covariate:

  survfit_oper <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + L,
                        data = trans_data[at_risk_oper == 1], cluster = ID,
                        timefix = timefix)
  survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + L,
                         data = trans_data, cluster = ID,
                         timefix = timefix)
  survfit_cens <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + L,
                        data = trans_data, cluster = ID,
                        timefix = timefix)
  survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0,
                       data = trans_data[at_risk_cov == 1], cluster = ID,
                       timefix = timefix)

  results[,2] <- c(survfit_oper$coefficients, survfit_death$coefficients[1], NA,
                   survfit_death$coefficients[2], survfit_cens$coefficients[1], NA,
                   survfit_cens$coefficients[2], survfit_cov$coefficients, NA)

  # Cox models only with operation:
  survfit_oper <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0,
                        data = trans_data[at_risk_oper == 1], cluster = ID,
                        timefix = timefix)
  survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A,
                         data = trans_data, cluster = ID,
                         timefix = timefix)
  survfit_cens <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A,
                        data = trans_data, cluster = ID,
                        timefix = timefix)
  survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + A,
                       data = trans_data[at_risk_cov == 1], cluster = ID,
                       timefix = timefix)

  results[,3] <- c(survfit_oper$coefficients, NA, survfit_death$coefficients, NA,
                   survfit_cens$coefficients, NA, survfit_cov$coefficients)

  return(results)

}


# Bootstrap for estimating model misspecification effect

sys_invest <- function(eff_L_A = 1, eff_L_D = 1, eff_A_D = -1, B = 200, N = 400) {

  # create results table
  mod_est <- matrix(nrow = 10*B, ncol = 3)

  # run simulations
  for(i in 1:B){
    #browser()
    sim <- sim_data_setting1(N = N, eta = rep(0.05,4), nu = rep(1.02, 4), 
                             beta_L_A = eff_L_A, beta_L_D = eff_L_D, beta_A_D = eff_A_D)
    mod_est[(i*10 - 9):(i*10),] <- fit_mod_set1(sim)
  }

  # create plots
  plotdata <- data.table(mod_est)
  names(plotdata) <- c("Fully Specified", "-A", "-L")

  plotdata[, "Est" := rep(c("Op L0", "Op L", "Death L0", "Death A", "Death L",
                            "Cens L0", "Cens A", "Cens L", "Cov L0", "Cov A"), B)]

  plotdata[, "True" := rep(c(1, eff_L_A, 1, eff_A_D, eff_L_D, 0, 0, 0, 1, -0.5), B)]

  b1 <- plotdata[,.("Bias True" = mean(`Fully Specified` - True)), by = Est]
  b2 <- plotdata[,.("Bias -A" = mean(`-A` - True)), by = Est]
  b3 <- plotdata[,.("Bias -L" = mean(`-L`- True)), by = Est]

  bias <- cbind(b1, b2[,2], b3[,2])

  pp1 <- ggplot(data = plotdata)+
    geom_histogram(aes(x = `Fully Specified`, y = ..density..), binwidth = 0.1,
                   fill = "blue", alpha = 0.5, col = "white")+
    facet_wrap(~Est, scales = "free") +
    geom_vline(aes(xintercept = True), color = "red", linetype = 2)

  pp2 <- ggplot(data = plotdata)+
    geom_histogram(aes(x = `-A`, y = ..density..), binwidth = 0.1,
                   fill = "blue", alpha = 0.5, col = "white")+
    facet_wrap(~Est, scales = "free") +
    geom_vline(aes(xintercept = True), color = "red", linetype = 2)

  pp3 <- ggplot(data = plotdata)+
    geom_histogram(aes(x = `-L`, y = ..density..), binwidth = 0.1,
                   fill = "blue", alpha = 0.5, col = "white")+
    facet_wrap(~Est, scales = "free") +
    geom_vline(aes(xintercept = True), color = "red", linetype = 2)

  return(list(Bias = bias, pp1 = pp1, pp2 = pp2, pp3 = pp3))

}


# Function for calculating the two probability fractions

prob_fracs <- function(data, tau = 5){
  
  data <- data[Time <= tau]
  
  # Individuals that experience T2D
  ID_T2D <- data[Delta == 3]$ID
  # Individuals that do not experience T2D
  ID_not_T2D <- unique(data[!ID %in% ID_T2D]$ID)
  
  # Number of individuals that experience T2D
  n1 <- ID_T2D |> length()
  # Number of individuals that do not experience T2D
  n0 <- ID_not_T2D  |> length()
  
  # Individuals that experience death with T2D 
  n3 <- ID_T2D  %in% data[Delta == 1]$ID |> sum()
  # Individuals that experience death without T2D 
  n2 <- ID_not_T2D  %in% data[Delta == 1]$ID |> sum()
  
  # Probability of experiencing T2D
  P_X1 <- n1 / N
  # Probability of not experiencing T2D
  P_X0 <- 1 - P_X1
  
  # Probability of dying with T2D
  P_X3 <- n3 / N
  # Probability of  of dying without T2D
  P_X2 <- n2 / N
  
  return(c("WO T2D" = P_X2 / P_X0,"W T2D" = P_X3 / P_X1))
}


# Function for finding the rolling average

roll_avg <- function(x) {
  n <- length(x)
  returnn <- numeric(n)
  for(i in 1:n) {
    returnn[i] <- mean(x[1:i])
  }
  returnn
}


