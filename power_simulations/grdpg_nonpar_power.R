
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
set.seed(1357)
MCs <- 500
epsilons <- c(0,.05,.1,.15,.2,.25,.3,.35,.4)
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
    vals[[i]][[j]] <- rep(0,MCs)
    k <- 1
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
                                          , Q = signs[[l]],numReps = 50)
      }
      
      cs <- sapply(get_matched,`[[`,3)
      Q <-  get_matched[[which.min(cs)]]$Q
      Xnew <- Xhat %*% Q
      vals[[i]][[j]][k] <-  nonpar.test(Xnew,Yhat,100)
      k <- k+1
    }
    #vals[[i]][[j]] <- vals[[i]][[j]]#/MCs
    j <- j + 1
  }
  names(vals[[i]]) <- ns
  
  i <- i + 1
}

names(vals) <- epsilons
save(vals,file = "MC_results_power_2.Rdata")

#code to create table in paper  
load("MC_results_power_2.Rdata") 
#vals2 <- vals
#load("simulation_results/5-2_GOOD_results/MC_results.Rdata") 

results <- matrix(0,9,5)
colnames(results) <- c(100,200,300,400,500)
rownames(results) <- names(vals)
results2 <- t(results)
rm(results)
results2[,"0"] <- sapply(vals$`0`,function(x){sum(ifelse(x < .1,1,0)) })
results2[,"0.05"] <- sapply(vals$`0`,function(x){sum(ifelse(x < .1,1,0)) })
results2[,"0.1"] <- sapply(vals$`0.1`,function(x){sum(ifelse(x < .1,1,0)) })
results2[,"0.15"] <- sapply(vals$`0.15`,function(x){sum(ifelse(x < .1,1,0)) })
results2[,"0.2"] <- sapply(vals$`0.2`,function(x){sum(ifelse(x < .1,1,0)) })
results2[,"0.25"] <- sapply(vals$`0.25`,function(x){sum(ifelse(x < .1,1,0)) })
results2[,"0.3"] <- sapply(vals$`0.3`,function(x){sum(ifelse(x < .1,1,0)) })
results2[,"0.35"] <- sapply(vals$`0.35`,function(x){sum(ifelse(x < .1,1,0)) })
results2[,"0.4"] <- sapply(vals$`0.4`,function(x){sum(ifelse(x < .1,1,0)) })


t(results2)
#make power curve
library(ggplot2)
library(reshape2)

gg_dat <- as.data.frame(t(results2/500))
gg_dat$eps <- colnames(results2)
gg_dat <- melt(gg_dat)
gg_dat$n <- as.integer(as.character(gg_dat$variable))
gg_dat$variable <- NULL
g <- ggplot(data=gg_dat,aes(x = n, y = value)) 
g+  theme_bw() +
  theme(plot.title = element_text(size = 10,hjust = 0.5))+
  ylab("Estimated Power") + 
  ggtitle("Power Curve for Different Values of Epsilon") +
  geom_line(aes(color = eps, group = eps),lwd=1,
            position=position_jitter(w=0.02, h=.012)) +
  geom_hline(yintercept = .05,linetype = "dashed")+
 # scale_x_continuous(limits = c(50,550)) +
  scale_y_continuous(limits = c(-0.02,1.02),breaks= c(0.00,0.05,0.25,0.5,.75,1.00)) 
  


load("MC_results_power_2.Rdata") 
#vals2 <- vals
#load("simulation_results/5-2_GOOD_results/MC_results.Rdata") 

results <- matrix(0,9,5)
colnames(results) <- c(100,200,300,400,500)
rownames(results) <- names(vals)
results2 <- t(results)
rm(results)
results2[,"0"] <- sapply(vals$`0`,mean)
results2[,"0.05"] <- sapply(vals$`0`,mean)
results2[,"0.1"] <- sapply(vals$`0.1`,mean)
results2[,"0.15"] <- sapply(vals$`0.15`,mean)
results2[,"0.2"] <- sapply(vals$`0.2`,mean)
results2[,"0.25"] <- sapply(vals$`0.25`,mean)
results2[,"0.3"] <- sapply(vals$`0.3`,mean)
results2[,"0.35"] <- sapply(vals$`0.35`,mean)
results2[,"0.4"] <- sapply(vals$`0.4`,mean)
results2

