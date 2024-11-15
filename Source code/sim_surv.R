sur_sim <- function(N,                      # Number of individuals
                    beta,                   # Effects
                    eta,                    # Shape parameters
                    nu,                     # Scale parameters
                    at_risk,                # Function defining the setting
                    term_deltas             # Terminal events
                    ){
  
  # Events
  x <- 1:ncol(beta)
  
  # Intensities
  phi_x <- function(x, k) {
    exp(L0 * beta[1,x] + k * beta[2,x] + A * beta[3,x] + L1[k+1] * beta[4,x])
  }
  
  lambda_x <- function(x, t, k, m) {
    at_risk(x - 1, k, m) * eta[x] * nu[x] * t ^ (nu[x] - 1) * phi_x(x, k)
  } 
  
  # Summed cumulative hazard
  sum_cum_haz <- function(u, t, k, m) {
    sum(sapply(x, function(x) {
      at_risk(x - 1, k, m) * eta[x] * phi_x(x, k) * ((t + u) ^ nu[x] - t ^ nu[x])
    }))}
  
  
  # Inverse summed cumulative hazard function
  inverse_sc_haz <- function(p, t, k, m, lower_bound = 10^-15, upper_bound = 100) {
    root_function <- function(u) sum_cum_haz(u, t, k, m) - p
    uniroot(root_function, lower = lower_bound, upper = upper_bound)$root
  }
  
  # Event probabilities
  probs <- function(t, k, m){
    probs <- sapply(x, function(x) lambda_x(x, t, k, m))
    summ <- sum(probs)
    probs / summ
  }
  
  # Simulation
  res <- data.frame(ID = numeric(), Time = numeric(), Delta = numeric(), L0 = numeric(), 
                L1 = numeric(), A = numeric())
  
  for(j in 1:N){
    # Draw
    L0 <- runif(1)
    A <- rbinom(1, 1, 0.5)
    L1 <- numeric()
    L1[1] <- runif(1)
    
    # Initialize
    T_k1 <- 0
    Delta <- -1
    
    # Iterate
    Ts <- c()
    Deltas <- c()
    # Indicator for number of events
    k <- 0
    # Indicator for operation
    m <- 0
    
    
    while(! Delta %in% term_deltas){
      # Update L1
      if(k > 0) L1[k + 1] <- rnorm(1, 0.05 * L0 + L1[k], sqrt(0.01))
      # Simulate time
      V <- runif(1)
      W_k1 <- inverse_sc_haz(-log(V), T_k1, k, m)
      T_k1 <- T_k1 + W_k1
      Ts[k + 1] <- T_k1
      # Simulate event
      Deltas[k + 1] <- sample(x, size = 1, prob = probs(T_k1, k, m)) - 1
      Delta <- Deltas[k + 1]
      # Update number of events and operation indicator
      k <- k + 1
      if(Delta == 3) m <- 1
    }
    
    jth_res <- data.frame(ID = rep(j, length(Ts)), 
                      Time = Ts, 
                      Delta = Deltas,
                      L0 = rep(L0,length(Ts)),
                      L1 = L1,
                      A = rep(A,length(Ts))) 
    res <- rbind(res, jth_res)
  }
  return(res)
}

