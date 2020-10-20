
#if(!require(devtools)) {
#  install.packages("devtools")
#  library(devtools)
#}

if(!require(nonparGraphTesting)) {
  install.packages("nonparGraphTesting_0.1.0.tar.gz", repos = NULL, type="source")
  library(nonparGraphTesting)
  
}
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

# functions to generate the graph given an arbitrary P matrix
# generateAdjacency simulates according to a P matrix
# generateMatrix creates the P matrix given b, the n dimensional assignment vector
# n, and bmatrix the B matrix


Rcpp::cppFunction("
  NumericMatrix generateAdjacencyMatrix(NumericMatrix pMatrix) {
  
    int n = pMatrix.cols();
    NumericMatrix A(n,n);
    for(int i = 0; i < n; i ++) {
      for (int j = i + 1; j < n; j++) {
        A(i,j) = (int)(rand()%100 < (pMatrix(i,j)* 100));
        A(j,i) = A(i,j);
      }
    }
    return A;
  }
")
set.seed(1111)
MCs <- 50
epsilons <- c(0,.05,.1,.15,.2,.25,.3)
ns <- c(100,200,300,400,500)#,1000)
vals <- list()
i <- 1
Q1 <- diag(1,3)
Q2 <- diag(-1,3)
Q3 <- diag(c(1,-1,-1),3)
Q4 <- diag(c(-1,1,1),3)
Q5 <- diag(c(-1,-1,1),3)
Q6 <- diag(c(-1,1,-1),3)
Q7 <- diag(c(1,1,-1),3)
Q8 <- diag(c(1,-1,1),3)
signs <- list(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8)
#n <- 100
for (eps in epsilons) {
  vals[[i]] <- list()
  print(paste0("epsilon = ",eps))
  B <- matrix(c(.3,.8,.8,.8,.3,.8,.8,.8,.3),3,3)
  B2 <- B + diag(eps,3)
  nus <- eigen(B)
  nus_true1 <- nus$vectors %*% diag(abs(nus$values)^(1/2),3,3)
  nus <- eigen(B2)
  nus_true2 <- nus$vectors %*% diag(abs(nus$values)^(1/2),3,3)
  Ipq <- diag(c(1,-1,-1),3,3)
  pis <- c(.35,.35,.3)
  j <- 1
  for (n in ns) {
    print(paste0("n = ",n))
    vals[[i]][[j]] <- 0
    for (mc in c(1:MCs)) {
      print(paste0("Run: ",mc," out of ", MCs, " for n = ",n, ", eps = ",eps))#Matching datasets")
      m <- n
      assignmentvector1 <- rmultinom(n,1,pis)
      assignmentvector2 <- rmultinom(m,1,pis)
      Xtrue <- t(assignmentvector1) %*% nus_true1
      Ytrue <- t(assignmentvector2) %*% nus_true2
      P1 <- Xtrue %*%Ipq %*% t(Xtrue)
      P2 <- Ytrue %*% Ipq %*% t(Ytrue)
      A <- generateAdjacencyMatrix(P1)
      C <- generateAdjacencyMatrix(P2)
      Xhat <- irlba(A,3)
      Xhat <- Xhat$u %*% diag(Xhat$d)^(1/2)
      #Xhat <- Xhat$vectors[,c(1,n-1,n)] %*% diag(abs(Xhat$values[c(1,n-1,n)])^(1/2))
      Yhat <- irlba(C,3)
      Yhat <- Yhat$u %*% diag(Yhat$d)^(1/2)
      #Yhat <- Yhat$vectors[,c(1,n-1,n)] %*% diag(abs(Yhat$values[c(1,n-1,n)])^(1/2))
      get_matched <- list()
      
      for (l in c(1:length(signs))) {
        get_matched[[l]] <- match_support(Xhat,Yhat
                                          ,lambda_init = .2
                                          ,alpha = .2
                                          , Q = signs[[l]],numReps = 0)
      }
      
      
      cs <- sapply(get_matched,`[[`,3)
      Q <-  get_matched[[which.min(cs)]]$Q
      Xnew <- Xhat %*% Q
      test <- nonpar.test(Xnew,Yhat,1000)
      vals[[i]][[j]] <-  sum(test$permutation_results) + vals[[i]][[j]]
    }
    vals[[i]][[j]] <- vals[[i]][[j]]/MCs
    j <- j + 1
  }
  names(vals[[i]]) <- ns
  
  i <- i + 1
}

names(vals) <- epsilons
#save(vals,file = "MC_results_9-29.Rdata")

#code to create table in paper  
#MCs <- 50
#epsilons <- epsilons
#ns <- c(100,200,300,400,500)#,1000)

results <- matrix(0,length(epsilons),length(ns))
colnames(results) <- ns
rownames(results) <- epsilons


results2 <- t(results)
# rm(results)
 for (eps in epsilons) {
   for(n in ns) {
     results2[as.character(n),as.character(eps)] <- vals[[as.character(eps)]][[as.character(n)]]
   }
 }

save(results2,file = "MC_results_9-30.Rdata")
# results2[,"0.1"] <- sapply(vals$`0.1`,`[[`,1)
# results2[,"0.15"] <- sapply(vals2$`0.15`,`[[`,1)
# results2[,"0.2"] <- sapply(vals2$`0.2`,`[[`,1)
# results2[,"0.25"] <- sapply(vals$`0.25`,`[[`,1)
# results2[,"0.3"] <- sapply(vals$`0.3`,`[[`,1)
# results2[,"0.35"] <- sapply(vals2$`0.35`,`[[`,1)
# round(t(results2/50),3)
