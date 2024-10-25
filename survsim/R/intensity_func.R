lambda <- function(x, t, k) {
  at_risk(x - 1, k) * eta[x] * nu[x] * t ^ (nu[x] - 1) * phi(x, k)
}
