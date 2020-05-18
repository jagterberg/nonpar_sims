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


#youngser's code to get elbow
getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main="") {
  ## Given a decreasingly sorted vector, return the given number of elbows
  ##
  ## Args:
  ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
  ##   n: the number of returned elbows.
  ##   threshold: either FALSE or a number. If threshold is a number, then all
  ##   the elements in d that are not larger than the threshold will be ignored.
  ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
  ##
  ## Return:
  ##   q: a vector of length n.
  ##
  ## Reference:
  ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
  ##   the scree plot via the use of profile likelihood", Computational
  ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006. 
  
  #  if (is.unsorted(-d))
  
  
  if (is.matrix(dat)) {
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
    d <- sort(dat,decreasing=TRUE)
  }
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev",main=main)
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b",main=main)
      points(q,dat[q],col=2,pch=19)
    }
  }
  
  return(q)
}


#load("data/TT-DSDS01216-glist114-raw-LCCTRUE.rda")
#dat <- list(A1s1,A2s1,A1s2,A2s2)
#save(dat,file = "data/dat.RData")
load("data/dat.RData")
A1s1 <- dat[[1]]
A2s1 <- dat[[2]]
A1s2 <- dat[[3]]
A2s2 <- dat[[4]]
rm(dat)
#diag(A1s1) <- rowSums(A1s1) / (nrow(A1s1)-1)
#diag(A2s2) <- rowSums(A2s2) / (nrow(A2s2)-1)
#diag(A1s2) <- rowSums(A1s2) / (nrow(A1s2)-1)
#diag(A2s1) <- rowSums(A2s1) / (nrow(A2s1)-1)


A1s1_eigen <- eigen(A1s1, symmetric = TRUE,only.values = TRUE)
A2s2_eigen <- eigen(A2s2, symmetric = TRUE,only.values=TRUE)
A1s2_eigen <- eigen(A1s2, symmetric = TRUE,only.values = TRUE)
A2s1_eigen <- eigen(A2s1, symmetric = TRUE,only.values = TRUE)



dA1s1 <- sort(abs(A1s1_eigen$values),decreasing = TRUE,index.return=TRUE)
dA2s2 <- sort(abs(A2s2_eigen$values),decreasing = TRUE,index.return=TRUE)
dA2s2 <- sort(abs(A2s2_eigen$values),decreasing = TRUE,index.return=TRUE)
dA1s2 <- sort(abs(A1s2_eigen$values),decreasing = TRUE,index.return=TRUE)


d1 <- getElbows(dA1s1$x,2)
d2 <- getElbows(dA2s2$x,2)
d3 <- getElbows(dA2s2$x,2)
d4 <- getElbows(dA1s2$x,2)
d <- max(d1,d2,d3,d4)
A1s1_svd <- irlba(A1s1,d)
A1s2_svd <- irlba(A1s2,d)
A2s1_svd <- irlba(A2s1,d)
A2s2_svd <- irlba(A2s2,d)
Xhat_A1s1 <- A1s1_svd$u[,c(1:d)] %*% diag(A1s1_svd$d[c(1:d)])^(1/2)
Xhat_A1s2 <- A1s2_svd$u[,c(1:d)] %*% diag(A1s2_svd$d[c(1:d)])^(1/2)
Xhat_A2s1 <- A2s1_svd$u[,c(1:d)] %*% diag(A2s1_svd$d[c(1:d)])^(1/2)
Xhat_A2s2 <- A2s2_svd$u[,c(1:d)] %*% diag(A2s2_svd$d[c(1:d)])^(1/2)
rm(A1s1,A1s2,A2s1,A2s2)

phat <- sum(ifelse(A1s1_eigen$values[dA1s1$ix[c(1:d)]] > 0,1,0))
qhat <- d - phat

#df.q <- expand.grid(d1=c(-1,1), d2=c(-1,1), d3=c(-1,1), d4=c(-1,1))
#QQ <- lapply(1:nrow(df.q), function(x) diag(df.q[x,]))

#out = NULL
#for (i in c(1:length(QQ))) {
#  out[[i]] <- match_support(X, Y, Q=QQ[[i]], numReps=50)
#}
#(costs = sapply(out, '[[', 3))
pq <- c(phat,qhat)
names(pq) <- c("phat","qhat")
results <- list()
results[[1]] <- pq

#A1s1 to A1s2
print("testing A1s1 to A1s2")
i <- 2
out <- match_support(Xhat_A1s1, Xhat_A1s2,  numReps=50)
pval <- nonpar.test(Xhat_A1s1 %*% out$Q ,Xhat_A1s2)
results[[i]] <- pval
names(results[[i]]) <- "A1s1 to A1s2"

#A1s1 to A2s1
print("testing A1s1 to A2s1")

i <- i+1
out <- match_support(Xhat_A1s1, Xhat_A2s1,  numReps=50)
pval <- nonpar.test(Xhat_A1s1  %*% out$Q ,Xhat_A2s1)
results[[i]] <- pval
names(results[[i]]) <- "A1s1 to A2s1"

#A1s1 to A2s2
print("testing A1s1 to A2s2")

i <- i+1
out <- match_support(Xhat_A1s1, Xhat_A2s2,  numReps=50)
pval <- nonpar.test(Xhat_A1s1  %*% out$Q ,Xhat_A2s2)
results[[i]] <- pval
names(results[[i]]) <- "A1s1 to A2s2"

#A1s2 to A2s1
print("testing A1s2 to A2s1")

i <- i+1
out <- match_support(Xhat_A1s2, Xhat_A2s1,  numReps=50)
pval <- nonpar.test(Xhat_A1s2  %*% out$Q ,Xhat_A2s1)
results[[i]] <- pval
names(results[[i]]) <- "A1s2 to A2s1"

#A1s2 to A2s2
print("testing A1s2 to A2s2")

i <- i+1
out <- match_support(Xhat_A1s2, Xhat_A2s2,  numReps=50)
pval <- nonpar.test(Xhat_A1s2  %*% out$Q ,Xhat_A2s2)
results[[i]] <- pval
names(results[[i]]) <- "A1s2 to A2s2"

#A2s1 to A2s2
print("testing A2s1 to A2s2 (last one)")
i <- i+1
out <- match_support(Xhat_A2s1, Xhat_A2s2,  numReps=50)
pval <- nonpar.test(Xhat_A2s1  %*% out$Q ,Xhat_A2s2)
results[[i]] <- pval
names(results[[i]]) <- "A2s1 to A2s2"
save(results,file ="real_data_results.Rdata")



