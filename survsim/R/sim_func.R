
# phi function
phi <- function(x, k) {
  exp(L0 * beta[,x][1] + k * beta[,x][2])
}


# intensity function
lambda <- function(x, t, k) {
  at_risk(x - 1, k) * eta[x] * nu[x] * t ^ (nu[x] - 1) * phi(x, k)
}


# sum of cumulative hazrds function
sum_cum_haz <- function(x, u, t, k) {
  sum(sapply(x, function(x) {
    at_risk(x - 1, k) * eta[x] * phi(x, k) * ((t + u) ^ nu[x] - t ^ nu[x])
  }))}



# inverse sum cumulative hazard function
inv_sum_cum_haz <- function(x, p, t, k, lower_bound = 0, upper_bound = 100) {

  root_function <- function(u) sum_cum_haz(x, u, t, k) - p

  uniroot(root_function, lower = lower_bound, upper = upper_bound)$root
}


# calculating probability of events
probs <- function(x, t, k){

  summ <- sum(sapply(x, function(x_i) lambda(x_i, t, k)))

  probs <- numeric(length(x))

  for(i in seq_along(x)){
    probs[i] <- lambda(x[i], t, k) / summ
  }

  return(probs)
}


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
      W_k1 <- inv_sum_cum_haz(x = x, p = -log(V), t = T_k1, K = k)
      T_k1 <- T_k1 + W_k1
      Ts[k + 1] <- T_k1
      Deltas[k + 1] <- sample(1:length(x), size = 1, prob = probs(x = x, t = T_k1, k = k)) - 1
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
