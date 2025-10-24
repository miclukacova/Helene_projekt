library(microbenchmark)

hely_metode <- function() {
  for (j in 1:num_events){
    event_times[,j] <- sapply(1:N, function(i) {
      lp_i_j <- lp_alive[[j]][i]
      basehaz_i_j <- basehaz_list[[j]][basehaz_list[[j]][, "time"]>T_k[i],]
      basehaz_i_j$haz <- cumsum(basehaz_i_j$dhazard*exp(lp_i_j))
      times_i_j <- basehaz_i_j$time[basehaz_i_j$haz<=U[i,j]]
      if(length(times_i_j)==0) {
        return(Inf) 
      } else {
        return(max(times_i_j))
      }
    })
  }
}

min_metode <- function() {
  # Calculate the Cox term
  cox_term <- lapply(cox_fits, 
                     function(model) exp(predict(model, newdata = sim_data, type="lp", reference = "zero")))
  
  # Calculate the cumulative intensity
  cumInt_list <- Map(function(base_vec, cox_vec) {
    outer(base_vec, cox_vec, "*")
  }, basehazz_list, cox_term)
  
  # We find the cumulative intensity of event j and the ith individual at the event time T_k
  cum_int_Tk <- matrix(ncol = num_events, nrow = num_alive)
  for(j in 1:num_events){
    cum_int_Tk[,j] <- sapply(1:num_alive, function(i) cum_int(cumInt_list[[j]][,i], jump_times[[j]], T_k[i]))
  }
  
  # Find the event times
  event_times <- matrix(ncol = num_events, nrow = num_alive)                  # matrix for event times
  for(j in 1:num_events){
    event_times[,j] <- sapply(1:num_alive, function(i) inv_cum_int(cumInt_list[[j]][,i], jump_times[[j]], V[i,j]))
  }
  
  # How many times can you experience the various events?
  for(j in 1:num_events){
    event_times[sim_data[, (2+j)] == n_event_max[j], j] <- Inf
  }
}


my_bench <- bench::mark(
      "hely" = hely_metode(),
      "mig" = min_metode(),
      check = FALSE
    )
)