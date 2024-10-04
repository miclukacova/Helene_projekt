#### Source code for simulation of survival data

# Libraries
library(gridExtra)
library(tidyverse)


















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
n <- 100
set.seed(98)
as <- runif(n)

# 4. Apply the recursive algorithm
ts <- numeric(n)
# i = 1:
ts[1] <- cum_hazard_inv( - log(as[1]))

for(i in 1 : (n - 1)){
  ts[i + 1] <- ts[i] + haz_tilde_inv( u = - log(as[i + 1]), t = ts[i])
}