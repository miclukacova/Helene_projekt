sum_cum_haz <- function(u, t, k) {
  sum(sapply(x, function(x) {
    at_risk(x - 1, k) * eta[x] * phi(x, k) * ((t + u) ^ nu[x] - t ^ nu[x])
  }))}
