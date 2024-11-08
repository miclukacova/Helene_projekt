sur_sim <- function(N,                      # Number of individuals
                    beta,                   # Baseline covariate effects
                    eta,                    # Shape parameters
                    nu,                     # Scale parameters
                    at_risk,                # Function defining the setting
                    term_deltas             # Terminal events
                    ){
  
  # Events
  x <- 1:ncol(beta)
  
  # Intensities
  phi_x <- function(x, k) {
    exp(L0 * beta[1,x] + k * beta[2,x] + A * beta[3,x])
  }
  
  
  lambda_x <- function(x, t, k) {
    at_risk(x - 1, k) * eta[x] * nu[x] * t ^ (nu[x] - 1) * phi_x(x, k)
  } 
  
  # Summed cumulative hazard
  sum_cum_haz <- function(u, t, k) {
    sum(sapply(x, function(x) {
      at_risk(x - 1, k) * eta[x] * phi_x(x, k) * ((t + u) ^ nu[x] - t ^ nu[x])
    }))}
  
  
  # Inverse summed cumulative hazard function
  inverse_sc_haz <- function(p, t, k, lower_bound = 0, upper_bound = 100) {
    
    root_function <- function(u) sum_cum_haz(u, t, k) - p
    
    uniroot(root_function, lower = lower_bound, upper = upper_bound)$root
  }
  
  # Event probabilities
  probs <- function(t, k){
    
    probs <- sapply(x, function(x) lambda_x(x, t, k))
    summ <- sum(probs)
    
    probs / summ
  }
  
  # Simulation
  res <- tibble(ID = numeric(), Time = numeric(), Delta = numeric(), L0 = numeric(), 
                L1 = numeric(), A = numeric())
  
  for(j in 1:N){
    #browser()
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
    k <- 0
    
    #browser()
    while(! Delta %in% term_deltas){
      if(k > 0) L1[k + 1] <- rnorm(1, 0.05 * L0 + L1[k], sqrt(0.01))
      V <- runif(1)
      W_k1 <- inverse_sc_haz(-log(V), T_k1, k)
      T_k1 <- T_k1 + W_k1
      Ts[k + 1] <- T_k1
      Deltas[k + 1] <- sample(1:length(x), size = 1, prob = probs(T_k1, k)) - 1
      Delta <- Deltas[k + 1]
      k <- k + 1
    }
    
    jth_res <- tibble(ID = rep(j, length(Ts)), 
                      Time = Ts, 
                      Delta = Deltas,
                      L0 = rep(L0,length(Ts)),
                      L1 = L1,
                      A = rep(A,length(Ts))) 
    res <- bind_rows(res, jth_res)
  }
  return(res)
}

