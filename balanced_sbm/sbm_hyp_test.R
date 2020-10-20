
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


run_simulation_sbm <- function(n=300,ntimes=100,seed=1234,eps=0,nMC=500) {
  results <- list()
  print(paste("initializing simulation for n =",n,"eps = ",eps))
  set.seed(seed) #1111 and #1112 is okay
  Q1 <- diag(1,2)
  Q2 <- diag(-1,2)
  Q3 <- diag(c(1,-1),2)
  Q4 <- diag(c(-1,1),2)
  signs <- list(Q1,Q2,Q3,Q4)#,Q5,Q6,Q7,Q8)
  signs2 <- list(as.matrix(1),as.matrix(-1))
  m <- n
  B <- matrix(c(.5,.8,.8,.8,.5,.8,.8,.8,.5),3,3)
  if(eps != 0) {
    B2 <- matrix(c(.5,.8,.8,.8,.5,.8,.8,.8,.5),3,3) + diag(eps,3)
  } else {
    B2 <- B
  }
  nus1 <- eigen(B)
  nus_true1 <- nus1$vectors %*% diag(abs(nus1$values)^(1/2),3,3)
  nus2 <- eigen(B2)
  nus_true2 <- nus2$vectors %*% diag(abs(nus2$values)^(1/2),3,3)
  Ipq <- diag(c(1,-1,-1),3,3)
  pis <- c(1/3,1/3,1/3)
  sigma <- 1/2
  #n <- 100
  
  print(paste("beginning simulation for n =",n,"eps = ",eps))
  for (i in c(1:ntimes)) {
    print(paste0("i = ",i," out of ",ntimes))
    assignmentvector1 <- rmultinom(n,1,pis)
    assignmentvector2 <- rmultinom(m,1,pis)
    Xtrue <-t(assignmentvector1) %*% nus_true1
    Ytrue <-  t(assignmentvector2) %*% nus_true2
    P1 <- Xtrue %*%Ipq %*% t(Xtrue)
    P2 <- Ytrue %*% Ipq %*% t(Ytrue)
    A <- generateAdjacencyMatrix(P1)
    C <- generateAdjacencyMatrix(P2)
    Xhat <- irlba(A,3)
    Xhat <- Xhat$u %*% diag(Xhat$d)^(1/2)
    Yhat <- irlba(C,3)
    Yhat <- Yhat$u %*% diag(Yhat$d)^(1/2)
    Xtilde <- irlba(P1,3)
    Xtilde <- Xtilde$u %*% diag(Xtilde$d)^(1/2)
    Ytilde <- irlba(P2,3)
    Ytilde <- Ytilde$u %*% diag(Ytilde$d)^(1/2)
    #find the alignment
    cs1 <- rep(0,length(signs))
    cs2 <- rep(0,length(signs))
    get_matched_1 <- list()
    get_matched_2 <- list()
    for (l in c(1:(length(signs)))) {
      get_matched_1[[l]] <- iterative_optimal_transport(Xhat,Yhat
                                                        #,lambda_init = .5
                                                        #,alpha = .5
                                                        #,lambda_final = .27
                                                        , Q = bdiag(1,signs[[l]]),numReps = 10
                                                        ,p=1,q=2)
      cs1[l] <- get_matched_1[[l]]$obj.value
      cs1[l] <- kernel.stat(Xhat%*% get_matched_1[[l]]$Q,Yhat)
      
      get_matched_2[[l]] <- iterative_optimal_transport(Xhat,Yhat
                                                        # ,lambda_init = .5
                                                        #,alpha = .5
                                                        # ,lambda_final = .27
                                                        , Q = bdiag(-1,signs[[l]]),numReps = 10
                                                        ,p=1,q=2)
      cs2[l] <- get_matched_2[[l]]$obj.value
      cs2[l] <- kernel.stat(Xhat%*% get_matched_2[[l]]$Q,Yhat)
    }
    
    minval1 <- cs1[which.min(cs1)]
    minval2 <- cs2[which.min(cs2)]
    
    
    W1 <- procrustes(Xhat,Xtilde,Pi = diag(1,n),p=1,q=2)
    W2 <- procrustes(Yhat,Ytilde,Pi = diag(1,n),p=1,q=2)
    
    Q1 <- iterative_optimal_transport(Xtilde,Ytilde,lambda=.01
                                      # ,lambda_init = .5
                                      #,alpha = .5
                                      # ,lambda_final = .27
                                      
                                      ,numReps =  50
                                      ,eps = .01,eps_OT = .01
                                      ,p=1,q=2)
    
    #Q1 <- Q1$Q
    
    Q_init <- W1 %*% Q1$Q %*% t(W2)
    
    
    get_matched_3 <- iterative_optimal_transport(Xhat,Yhat,lambda=.01
                                                 # ,lambda_init = .5
                                                 #,alpha = .5
                                                 # ,lambda_final = .27
                                                 , Q = Q_init
                                                 ,numReps =  50
                                                 ,eps = .01,eps_OT = .01
                                                 ,p=1,q=2)
    
    minval3 <- kernel.stat(Xhat%*% get_matched_3$Q,Yhat)

    
    if( minval1 < minval2 & minval1 < minval3) { 
      final_Q <- get_matched_1[[which.min(cs1)]]$Q
      
     
    } else if ( minval3 < minval1){
      final_Q <- get_matched_3$Q#[[which.min(cs1)]]$Q
    } else {
      final_Q <- get_matched_2[[which.min(cs2)]]$Q
    }
    
    Xnew <- Xhat%*% final_Q
    
    results[[i]] <- nonpar.test(Xnew,Yhat,nsims = nMC)
  }
  
  
  
  alpha <- .05
  for (i in c(1:length(results))) {
    results[[i]]$critical_value <- sort(results[[i]]$permutation_results)[floor((1-alpha)*length(results[[i]]$permutation_results))]
    results[[i]]$`estimated p-value` <- sum(results[[i]]$permutation_results > results[[i]]$`test statistic`)/length(results[[i]]$permutation_results)
  }
  mainvals <- paste0("n = ",n,", eps = ",eps)
  toReturn <-list(mainvals,results)
  return(toReturn)
  
}














