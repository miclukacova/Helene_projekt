library(foreach)
library(doParallel)
library(itertools)
library(randomForestSRC)

num_splits <-4
cl <- makeCluster(num_splits)
registerDoParallel(cl)

predictions <-
  foreach(d=isplitRows(sim_data, chunks=4), 
          .packages=c("randomForestSRC")) %dopar% {
            predict(RF_fit, newdata=d)
          }

stopCluster(cl)

# Combine results
chf_mat <- combined <- do.call(abind::abind, c(predictions, along = 1))
times <-  c(0,pred_list[[1]]$time.interest)                                   # Time points
