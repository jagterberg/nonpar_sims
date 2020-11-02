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



# Simulate two SBMs and test
run_simulation_dcsbm <- function(n=300,ntimes=100,seed=2020,nMC=500,betaparams =c(.2,.8,1,1)) {
  print(paste("beginning simulation for n =",n))
  alpha <- .05
  results <- list()
  set.seed(seed) #1111 and #1112 is okay
  Q1 <- diag(1,2)
  Q2 <- diag(-1,2)
  Q3 <- diag(c(1,-1),2)
  Q4 <- diag(c(-1,1),2)
  signs <- list(Q1,Q2,Q3,Q4)#,Q5,Q6,Q7,Q8)
  signs2 <- list(as.matrix(1),as.matrix(-1))
  m <- n
  B <- matrix(c(.5,.8,.8,.8,.5,.8,.8,.8,.5),3,3)
  nus <- eigen(B)
  nus_true <- nus$vectors %*% diag(abs(nus$values)^(1/2),3,3)
  Ipq <- diag(c(1,-1,-1),3,3)
  pis <- c(1/3,1/3,1/3)
  sigma <- 1/2
 
  #n <- 100
  for (i in c(1:ntimes)) {
    print(paste0("i = ",i," out of ",ntimes))
    assignmentvector1 <- rmultinom(n,1,pis)
    assignmentvector2 <- rmultinom(m,1,pis)
    betas <- betaparams[1]*rbeta(n,betaparams[3],betaparams[4]) + betaparams[2]
    Xtrue <-t(assignmentvector1) %*% nus_true
    Ytrue <-  betas* t(assignmentvector2) %*% nus_true
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
    cs3 <- rep(0,length(signs))
    cs4 <- rep(0,length(signs))
    cs5 <- rep(0,length(signs))
    cs6 <- rep(0,length(signs))
    get_matched_1 <- list()
    get_matched_2 <- list()
    get_matched_5 <- list()
    get_matched_6 <- list()
    gm <- list()
    gm2 <- list()
    for (l in c(1:(length(signs)))) {
      get_matched_1[[l]] <- iterative_optimal_transport(Xhat,Yhat,lambda=.001
                                                        #,lambda_init = .5
                                                        #,alpha = .5
                                                        #,lambda_final = .27
                                                        , Q = bdiag(1,signs[[l]]),numReps = 10
                                                        ,p=1,q=2)
      #cs1[l] <- get_matched_1[[l]]$obj.value
      cs1[l] <- kernel.stat(Xhat%*% get_matched_1[[l]]$Q,Yhat)
      
      get_matched_2[[l]] <- iterative_optimal_transport(Xhat,Yhat,lambda=.001
                                                        # ,lambda_init = .5
                                                        #,alpha = .5
                                                        # ,lambda_final = .27
                                                        , Q = bdiag(-1,signs[[l]]),numReps = 10
                                                        ,p=1,q=2)
      #cs2[l] <- get_matched_2[[l]]$obj.value
      cs2[l] <- kernel.stat(Xhat%*% get_matched_2[[l]]$Q,Yhat)
      gm[[l]] <- iterative_optimal_transport(Xtilde,Ytilde,lambda=.001
                                             # ,lambda_init = .5
                                             #,alpha = .5
                                             # ,lambda_final = .27
                                             , Q = bdiag(-1,signs[[l]]),numReps = 10
                                             ,p=1,q=2)
      cs3[l] <- kernel.stat(Xtilde%*% gm[[l]]$Q,Ytilde)
      
      gm2[[l]] <- iterative_optimal_transport(Xtilde,Ytilde,lambda=.001
                                              # ,lambda_init = .5
                                              #,alpha = .5
                                              # ,lambda_final = .27
                                              , Q = bdiag(1,signs[[l]]),numReps = 10
                                              ,p=3)#,q=2)
      cs4[l] <- kernel.stat(Xtilde%*% gm2[[l]]$Q,Ytilde)
      
      get_matched_5[[l]] <- iterative_optimal_transport(Xhat,Yhat,lambda=.001
                                                        #,lambda_init = .5
                                                        #,alpha = .5
                                                        #,lambda_final = .27
                                                        , Q = bdiag(1,signs[[l]]),numReps = 10
                                                        ,p=3)#,q=2)
      #cs1[l] <- get_matched_1[[l]]$obj.value
      get_matched_5[[l]]$Q <- procrustes(X=diag(1,3),Y=get_matched_5[[l]]$Q,Pi = diag(1,3),p=1,q=2)
      
      cs5[l] <- kernel.stat(Xhat%*% get_matched_5[[l]]$Q,Yhat)
      
      get_matched_6[[l]] <- iterative_optimal_transport(Xhat,Yhat,lambda=.001
                                                        # ,lambda_init = .5
                                                        #,alpha = .5
                                                        # ,lambda_final = .27
                                                        , Q = bdiag(-1,signs[[l]]),numReps = 10
                                                        ,p=3)#1,q=2)
      #cs2[l] <- get_matched_2[[l]]$obj.value
      get_matched_6[[l]]$Q <- procrustes(X=diag(1,3),Y=get_matched_6[[l]]$Q,Pi = diag(1,3),p=1,q=2)
      cs6[l] <- kernel.stat(Xhat%*% get_matched_6[[l]]$Q,Yhat)
      
      gm[[l]] <- iterative_optimal_transport(Xtilde,Ytilde,lambda=.001
                                             # ,lambda_init = .5
                                             #,alpha = .5
                                             # ,lambda_final = .27
                                             , Q = bdiag(-1,signs[[l]]),numReps = 10
                                             ,p=1,q=2)
      cs3[l] <- kernel.stat(Xtilde%*% gm[[l]]$Q,Ytilde)
      
      
    }
    
    minval1 <- cs1[which.min(cs1)]
    minval2 <- cs2[which.min(cs2)]
    minval3 <- cs3[which.min(cs3)]
    minval4 <- cs4[which.min(cs4)]
    minval5 <- cs5[which.min(cs5)]
    minval6 <- cs6[which.min(cs6)]
    
    
    #get_q_init <- cs3[which.min(cs3)]
    W1 <- procrustes(Xhat,Xtilde,Pi = diag(1,n),p=1,q=2)
    W2 <- procrustes(Yhat,Ytilde,Pi = diag(1,n),p=1,q=2)
    
    gms <- c(cs3,cs4)
    
    mm <- which.min(gms)
    
    if (mm > length(signs)) {
      Q0 <- gm[[mm- length(signs)]]$Q
    } else {
      Q0 <- gm2[[mm]]$Q
    }
    
    Q1 <- iterative_optimal_transport(Xtilde,Ytilde,lambda=.001
                                      # ,lambda_init = .5
                                      #,alpha = .5
                                      # ,lambda_final = .27
                                      ,Q = Q0
                                      ,numReps =  10
                                      ,eps = .001,eps_OT = .001
                                      ,p=1,q=2)
    
    
    
    #Q1 <- Q1$Q
    #Q_indef <- get_initialization(Xtilde,Ytilde,p=1,q=2)
    W_X <- procrustes(Xhat, Xtilde,Pi=diag(1,n))
    W_Y <- procrustes(Yhat, Ytilde,Pi=diag(1,n))
    #W_tilde <- procrustes(Xtilde,Ytilde,Pi=diag(1,n))
    
    
    Q_init1 <- W1 %*% Q1$Q %*% t(W2)
    Q_init2 <- W_X %*% Q1$Q  %*% t(W_Y)
    
    get_matched_3 <- iterative_optimal_transport(Xhat,Yhat,lambda=.001
                                                 # ,lambda_init = .5
                                                 #,alpha = .5
                                                 # ,lambda_final = .27
                                                 , Q = Q_init1
                                                 ,numReps =  10
                                                 ,eps = .001,eps_OT = .001
                                                 ,p=1,q=2)
    get_matched_4 <- iterative_optimal_transport(Xhat,Yhat,lambda=.001
                                                 # ,lambda_init = .5
                                                 #,alpha = .5
                                                 # ,lambda_final = .27
                                                 , Q = Q_init2
                                                 ,numReps =  10
                                                 ,eps = .001,eps_OT = .001
                                                 ,p=3)#,q=2)
    
    get_matched_4$Q <- procrustes(X=diag(1,3),Y=get_matched_4$Q,Pi = diag(1,3),p=1,q=2)
    minval3 <- kernel.stat(Xhat%*% get_matched_3$Q,Yhat)
    minval4 <- kernel.stat(Xhat%*%Q_init1,Yhat)
    minval7 <- kernel.stat(Xhat%*%get_matched_4$Q,Yhat)
    
    w.min <- which.min(c(minval1,minval2,minval3,minval4,minval5,minval6,minval7))
    
    
    if( w.min==1) { 
      final_Q <- get_matched_1[[which.min(cs1)]]$Q
    } else if( w.min  ==2) {
      final_Q <- get_matched_2[[which.min(cs2)]]$Q
    } else if ( w.min  ==3) {
      final_Q <- get_matched_3$Q#[[which.min(cs1)]]$Q
    } else if ( w.min  ==4) {
      final_Q <- Q_init1
    } else if (  w.min  ==5){
      final_Q <- get_matched_5[[which.min(cs5)]]$Q
    }else if (  w.min  ==6){
      final_Q <- get_matched_6[[which.min(cs6)]]$Q
    } else {
      final_Q <- get_matched_4$Q
    }
    
    
    
    
    
    Xnew <- Xhat%*% final_Q
    
    results[[i]] <- nonpar.test(Xnew,Yhat,nsims = nMC)
  }
  
  alpha <- .05
  for (i in c(1:length(results))) {
    results[[i]]$critical_value <- sort(results[[i]]$permutation_results)[floor((1-alpha)*length(results[[i]]$permutation_results))]
    results[[i]]$`estimated p-value` <- sum(results[[i]]$permutation_results > results[[i]]$`test statistic`)/length(results[[i]]$permutation_results)
  }
  mainvals <- paste0("n = ",n)
  toReturn <-list(mainvals,results)
  print(paste("finished simulation for n =",n))
  return(toReturn)
  
}
