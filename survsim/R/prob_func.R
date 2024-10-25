probs <- function(t, k){

  summ <- sum(sapply(x, function(x_i) lambda(x_i, t, k)))

  probs <- numeric(length(x))

  for(i in seq_along(x)){
    probs[i] <- lambda(x[i], t, k) / summ
  }

  return(probs)
}
