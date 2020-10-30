if(!require(foreach)) {
  install.packages("foreach")
  library(foreach)
}
if(!require(parallel)) {
  install.packages("parallel")
  library(parallel)
}
if(!require(doParallel)) {
  install.packages('doParallel')
  library(doParallel)
}
install.packages("nonparGraphTesting_0.1.0.tar.gz", repos = NULL, type="source")
library(nonparGraphTesting)
if (!require(irlba)) {
  install.packages("irlba")
  library(irlba)
}
if(!require(igraph)) {
  install.packages("igraph")
  library(igraph)
}
if(!require(Rcpp)) {
  install.packages("Rcpp")
  library(Rcpp)
}
if(!require(Matrix)) {
  install.packages("Matrix")
  library(Matrix)
}

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#epsilons <- c(0,.1,.2)
ns <- c(100,200,300,400,500,600,700)
bps <- list(
  c(0,1,1,1),c(.1,.9,1,1), c(.2,.8,1,1), c(.3,.7,1,1)
)
print(paste0("packages loaded, running dcSBM simulation"))#,cl))


results_dcsbm <- 
  foreach(n=ns,.packages=c('nonparGraphTesting','irlba','igraph','Rcpp','Matrix')
          ,.noexport = "generateAdjacencyMatrix") %:% 
  foreach(bs = c(1,2,3,4),.packages=c('nonparGraphTesting','irlba','igraph','Rcpp','Matrix')
           ,.noexport = "generateAdjacencyMatrix" )  %dopar% {
             source("./balanced_vs_dcsbm/sbm_vs_dcsbm.R")
             #print(paste("eps = ",eps,", n = ",n))
             run_simulation_dcsbm(ntimes = 100,n=n,nMC = 500,betaparams = bps[[bs]])
    }
save(results_dcsbm,file = "dcsbm_results_10-26.Rdata") 
stopCluster(cl)
print("finished.")

