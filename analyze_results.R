load("sbm_results_10-14.Rdata")

#extract estimated p-values
k <- 1
vals <- list()
for (i in 1:length(results_sbm)) {
  vals[[i]] <- matrix(0,3,500)
  for (j in 1:length(results_sbm[[i]])) {
    
    for (k in c(1:length(results_sbm[[i]][[j]])))
      vals[[i]][j,k] <- results_sbm[[i]][[j]][[k]]$`estimated p-value`
      
  }
}

power_estimate_300 <- 1 - rowSums(ifelse(vals[[1]] > .05,1,0))/500
power_estimate_600 <- 1 - rowSums(ifelse(vals[[2]] > .05,1,0))/500
power_estimate_900 <- 1 - rowSums(ifelse(vals[[3]] > .05,1,0))/500

powerv <- rbind(power_estimate_300,power_estimate_600,power_estimate_900)
rownames(powerv) <- c("300","600","900")
colnames(powerv) <- c("0",".1",".2")






