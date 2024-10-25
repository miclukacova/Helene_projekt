#' @export
sur_sim <- function(N,                      # Number of individuals
                    beta,                   # Baseline covariate effects
                    eta,                    # Shape parameters
                    nu,                     # Scale parameters
                    at_risk,                # Function defining the setting
                    term_deltas             # Terminal events
){

  x <- 1:ncol(beta)

  # Table for results
  res <- tibble(ID = numeric(), L0 = numeric(), Time = numeric(), Delta = numeric())

  for(j in 1:N){

    # Draw baseline covariates
    L0 <- runif(1)

    # Initialize
    T_k1 <- 0
    Delta <- -1

    # Iterate
    Ts <- c()
    Deltas <- c()
    k <- 0

    # While we are not in a terminal event
    while(! Delta %in% term_deltas){
      V <- runif(1)
      W_k1 <- inv_sum_cum_haz(-log(V), T_k1, k)
      T_k1 <- T_k1 + W_k1
      Ts[k + 1] <- T_k1
      Deltas[k + 1] <- sample(1:length(x), size = 1, prob = probs(T_k1, k)) - 1
      Delta <- Deltas[k + 1]
      k <- k + 1
    }

    jth_res <- tibble(ID = rep(j, length(Ts)),
                      L0 = rep(L0,length(Ts)),
                      Time = Ts,
                      Delta = Deltas)
    res <- bind_rows(res, jth_res)
  }
  return(res)
}
