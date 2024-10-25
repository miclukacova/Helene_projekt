inv_sum_cum_haz <- function(p, t, k, lower_bound = 0, upper_bound = 100) {

  root_function <- function(u) sum_cum_haz(t, k) - p

  uniroot(root_function, lower = lower_bound, upper = upper_bound)$root
}
