
# No treatment, L1 or k effect and no additional covariate process

simple1 <- function(N,                      # Number of individuals
                    beta,                   # Effect of L0 on operation, censoring, death
                    eta,                    # Shape parameters
                    nu                      # Scale parameters
                    ){
  
  at_risk <- function(x, k, m) {
    # If you have not died yet or been censored yet, you are at risk for dying or being censored
    if(x == 1 | x == 2) return(1)
    # You are only at risk for an operation if you have not had an operation yet
    else if(x == 0) return(as.numeric(k == 0 | (k == 1 & m == 1)))
    # You are never at risk for a change in the covariate process
    else return(0)
  }
  
  eta <- c(eta, 0)
  nu <- c(nu, 0)
  term_deltas <- c(1, 2)
  
  beta <- rbind(c(beta,0), rep(0,4), rep(0,4), rep(0,4))
  results <- sur_sim(N, beta, eta, nu, at_risk, term_deltas)
  
  return(results[,1:4])
}
  