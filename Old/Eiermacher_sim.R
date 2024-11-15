# Simulation based on Eiermacher paper
# 1. The Weibull hazard function as a function of total time

hazard <- function(t, lambda = 2, v = 0.5) lambda * v * t ^ (v - 1) 
cum_hazard <- function(t, lambda = 2, v = 0.5) lambda * t ^ (v)
cum_hazard_inv <- function(t, lambda = 2, v = 0.5) (t / lambda) ^ (1 / v)

# 2. Derive the cumulative hazard function for interevent time

haz_tilde <- function(u, t, lambda = 2, v = 0.5) {
  lambda * ((t + u) ^ v - t ^ v)
}

haz_tilde_inv <- function(u, t, lambda = 2, v = 0.5) {
  ((u + lambda * t ^ v) / lambda) ^ (1 / v) - t
}

# 3. Simulate independent random numbers a_i following a uniform distribution 

# Number of individuals
N <- 1000
times <- vector("list", length = N)
cens_times <- numeric(N)

# Number of events
n <- 100

set.seed(98)
for(j in 1:N){
  
  as <- runif(n)
  ts <- c()
  C <- rweibull(1, 35, 10)
  
  # 4. Apply the recursive algorithm
  i <- 1
  t <- cum_hazard_inv( - log(as[i]))
  while(t <= C){
    ts <- c(ts, t)
    i <- i + 1
    t <- t + haz_tilde_inv(u = - log(as[i]), t = t)
  }
  
  times[[j]] <- c(ts[ts < C])
  cens_times[j] <- C 
}

# Output 
id <- c()
time <- c()

for(i in 1:N){
  id <- c(id, rep(i, length(times[[i]])))
  time <- c(time, times[[i]])
}

tibble(Id = id, Time = time)


