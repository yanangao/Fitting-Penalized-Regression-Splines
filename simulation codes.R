
rm(list=ls())
# ============================Functions============================
# 1. True function f: customized???
# Input: x
# Return: y
true_f1 <- function(x) {
  # Define your own true functon
  y <- 1.5*dnorm((x-0.35)/0.15) - dnorm((x-0.8)/0.04)
  return(y)
}
true_f2 <- function(x,j) {
  temp <- 2^((9-4*j)/5)
  y <- sqrt(x*(1-x))*sin( (2*pi*(1+temp))/(x+temp) )
  return(y)
}
# 2. Sample X
# Input: sample size
# Return: a vector of sampled x, size = n
sample_x <- function(n,j) {
  shape1 <- (j+4)/5
  shape2 <- (11-j)/5
  p <- runif(n,0,1)
  x <- qbeta(p,shape1,shape2,lower.tail=TRUE)
}
# 3. Construct Smoothing Matrix
# Input: sampled x, lambda, spline degree p, # of knots nk
# Return: the smmoothing matrix H = (t(X)X+lambda*D)^{-1}*t(X)
H_lambda <- function(x,lambda,p,nk) {
  knots <- seq(min(x),max(x),(max(x)-min(x))/(1+nk))[c(-1,-(nk+2))]
  n_row <- length(x)
  n_col <- p+1+length(knots)
  # Initialize design matrix
  X <- matrix(1,nrow=n_row,ncol=n_col)
  for(k in seq(2,p+1)) {
    X[,k] <- x^(k-1)
  }
  for(k in seq(1,nk)) {
    tk <- knots[k]
    index <- ((x-tk) >= 0)*1
    X[,k+(p+1)] <- (index*(x-tk))^p
  }
  
  # temp matrix stores t(X)*X
  temp <- t(X)%*%X
  for(k in seq(1,nk,1)) {
    pos <- k+(p+1)
    temp[pos,pos] <- temp[pos,pos] + lambda
  }
  
  smoothing_matrix <- X%*%solve(temp,t(X))
  
  return(smoothing_matrix)
}

# 4. CV score
# Input: smoothing matrix S, observerd y
# Return: CV score
CV_score <- function(y,S) {
  S_diag <- diag(S)
  y_hat <- S%*%y
  CV <- mean(((y-y_hat)/(1-S_diag))^2)
  return(CV)
}

# 5. GCV score
# Inpute: smoothing matrix S, observed y
# Return: GCV score
GCV_score <- function(y,S) {
  S_diag <- diag(S)
  y_hat <- S%*%y
  n <- length(y)
  GCV <- mean((y - y_hat)^2)/(1-sum(S_diag)/n)^2
  return(GCV)
}

# 6. AICc score
# Input: smoothing matrix S, observered y
# Return: corrected AIC score
AICc_score <- function(y,S) {
  trace_S <- sum(diag(S))
  y_hat <- S%*%y
  n <- length(y)
  AICc <- log(sum((y-y_hat)^2)) + 2*(trace_S+1)/(n-trace_S-2)
  return(AICc)
}

# 7. Risk score
# Input: smoothing matrix S, observed y, sigma_sq_hat pre-chosen by cv
# Return: Risk score
# The estimation of sigma^2 follows from the reference paper on Page 141
Risk_score <- function(y,S,sigma_sq_hat){
  
  y_hat <- S%*%y
  n <- length(y)
  trace_S <- sum(diag(S))
  error_l2_norm <- sum((y-y_hat)^2)
  risk <- error_l2_norm - sigma_sq_hat*(n-2*trace_S)
  return(risk)
}

# 8. Get the optimal lambda
# Input: observed x, observed y, score function, spline degree p, 
#        # of knots nk, score function name score_fun
# Return: the optimal lambda
# By default, lambda is taken in 10^(-seq(10,-2,-0.05))
Optimal_Lambda <- function(x,y,p,nk,score_fun) {
  
  # Prepare scores 
  prepare_scores <- function(x,y,p,nk,lambda_vec,score_fun) {
    scoreFUN <- match.fun(score_fun)
    # Instead of using H_lambda function
    # Inside this function we construct design matrix X once
    # Use it to compute all score values
    knots <- seq(min(x),max(x),(max(x)-min(x))/(1+nk))[c(-1,-(nk+2))]
    n_row <- length(x)
    n_col <- p+1+length(knots)
    # Initialize design matrix
    X <- matrix(1,nrow=n_row,ncol=n_col)
    for(k in seq(2,p+1)) {
      X[,k] <- x^(k-1)
    }
    for(k in seq(1,nk)) {
      tk <- knots[k]
      index <- ((x-tk) >= 0)*1
      X[,k+(p+1)] <- (index*(x-tk))^p
    }
    
    # temp matrix stores t(X)*X
    temp <- t(X)%*%X
    
    # Continue the work HERE ......
    score_vec <- c()
    
    # Compute score for each lambda
    for(lambda in lambda_vec) {
      for(k in seq(1,nk,1)) {
        pos <- k+(p+1)
        temp[pos,pos] <- temp[pos,pos] + lambda
      }
      
      S <- X%*%solve(temp,t(X))
      
      score_vec <- c(score_vec, scoreFUN(y,S))
    }
    return(score_vec)
  }
  
  # Values that lambda can take
  lambda_vec <- 10^(-seq(10,-2,-0.1))
  n <- length(y)
  
  
  if(score_fun == 'Risk_score') {
    
    score_temp <- prepare_scores(x,y,p,nk,lambda_vec,'CV_score')
    optimal_lambda_temp <- lambda_vec[which(score_temp == min(score_temp))]
    
    S_cv <- H_lambda(x,optimal_lambda_temp,p,nk)
    y_hat_cv <- S_cv%*%y
    trace_S_cv <- sum(diag(S_cv))
    
    sigma_sq_hat <- sum((y - y_hat_cv)^2)/(n - trace_S_cv)
    
    scoreFUN <- function(y,S) { Risk_score(y,S,sigma_sq_hat) }
  } else {
    scoreFUN <- match.fun(score_fun)  
  }
  
  
  score_vec <- prepare_scores(x,y,p,nk,lambda_vec,scoreFUN)
  
  optimal_index <- which(score_vec == min(score_vec))
  
  return(lambda_vec[optimal_index])
}
# ============================SETUP============================
# Order of spline
p <- 3
# Number of knots
nk <- 30
# Model assumption: y = f + e, e iid mean 0, var sigma^2
# Simluation 
T <- 200 # 200
# Cases
J <- 6 # 6
# Sample size
n <- 200 # 200
# Wilcoxon Rank significance level
alpha <- 0.05/4
# ==========================Noise Level========================
# Save Plot
jpeg('noise-level.jpg', width=4.25, height=3.25, units="in", res=1000, pointsize=4)
par(mfrow=c(3,4))


for(j in seq(1,J)) {
  
  Lambda_matrix <- matrix(0,nrow=4,ncol=T)
  rownames(Lambda_matrix) <- c('CV','GCV','AIC','Risk')
  r_matrix <- matrix(0,nrow=4,ncol=T)
  rownames(r_matrix) <- c('CV','GCV','AIC','Risk')
  XY_matrix <- matrix(0,nrow=n,ncol=2*T)
  
  
  for(t in seq(1,T)) {
    set.seed(14067+j*t)
    # Normal error with mean 0 and var 1
    e <- rnorm(n, 0, 1)
    # sigma = 0.1
    e <- (0.02+0.04*(j-1)^2)*e
    # Sample x 
    x <- seq((1-0.5)/n,(n-0.5)/n,1/n)
    # Construct observerd y
    f <- true_f1(x)
    y <- f + e
    
    XY_matrix[,2*t-1] <- x
    XY_matrix[,2*t] <- y
    
    err_smallest <- unlist(optimize(function(z) { sum((f-H_lambda(x,z,p,nk)%*%y)^2) } ,c(0,100), tol=.Machine$double.eps))[2]
    
    opt_cv <- Optimal_Lambda(x,y,p,nk,'CV_score')
    opt_gcv <- Optimal_Lambda(x,y,p,nk,'GCV_score')
    opt_aicc <- Optimal_Lambda(x,y,p,nk,'AICc_score')
    opt_risk <- Optimal_Lambda(x,y,p,nk,'Risk_score')
    
    Lambda_matrix[1,t] <- opt_cv
    Lambda_matrix[2,t] <- opt_gcv
    Lambda_matrix[3,t] <- opt_aicc
    Lambda_matrix[4,t] <- opt_risk
    
    r_matrix[1,t] <- sum((f - H_lambda(x,opt_cv,p,nk)%*%y)^2)/err_smallest
    r_matrix[2,t] <- sum((f - H_lambda(x,opt_gcv,p,nk)%*%y)^2)/err_smallest
    r_matrix[3,t] <- sum((f - H_lambda(x,opt_aicc,p,nk)%*%y)^2)/err_smallest
    r_matrix[4,t] <- sum((f - H_lambda(x,opt_risk,p,nk)%*%y)^2)/err_smallest
    
  }
  
  r_matrix <- log(r_matrix)
  
  wilcox_rank <- vector(mode='list', length=4)
  names(wilcox_rank) <- c('CV', 'GCV', 'AICc', 'Risk')
  
  med_vec <- apply(r_matrix,1,median)
  
  rank_vec <- order(med_vec)
  
  test <- c()
  for(k in seq(1,3,1)) {
    r1 <- r_matrix[rank_vec[k],]
    r2 <- r_matrix[rank_vec[k+1],]
    ind <- (1-((r1-r2)==0)*1)==1
    if(sum(ind*1)==0) {
      test <- c(test,0)
    } else{
      test <- c(test,(wilcox.test(r1[ind],r2[ind],paired=T,alternative='less')$p.value < alpha)*1)
    }
  }
  
  if (sum(test) == 0) {
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- 2.5
    }
  } else if (sum(test) == 1) {
    pos <- which(test == 1)
    temp1 <- (1+pos)/2
    temp2 <- (pos+5)/2
    for(k in seq(1,pos)) {
      wilcox_rank[[rank_vec[k]]] <- temp1
    }
    for(k in seq(pos+1,4)) {
      wilcox_rank[[rank_vec[k]]] <- temp2
    }
  } else if (sum(test) == 2) {
    pos <- which(test == 0)
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- k
    }
    wilcox_rank[[rank_vec[pos]]] <- pos+0.5
    wilcox_rank[[rank_vec[pos+1]]] <- pos+0.5
  } else {
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- k
    }
  }
  
  print(unlist(wilcox_rank))
  
  # Plot true function
  plot(x,f,type='l',col='blue',main=paste('j=',j),xlab='',ylab='',ylim=c(-0.5,1.5))
  # Plot a typical simulation x,y
  points(XY_matrix[,1], XY_matrix[,2],col='red')
  
  boxplot_names <- c(wilcox_rank$CV,wilcox_rank$GCV,wilcox_rank$AICc,wilcox_rank$Risk)
  boxplot(r_matrix[1,],r_matrix[2,],r_matrix[3,],r_matrix[4,],names=boxplot_names,xlab='',ylab='',ylim=c(0,1.5))
  
  save(Lambda_matrix,XY_matrix,r_matrix,file=paste('noise-level-J',j,'.RData',sep=''))
  
}
dev.off()

# ==========================Design Density========================
jpeg('design-density.jpg', width=4.25, height=3.25, units="in", res=1000, pointsize=4)
par(mfrow=c(3,4))

for(j in seq(1,J)) {
  
  Lambda_matrix <- matrix(0,nrow=4,ncol=T)
  rownames(Lambda_matrix) <- c('CV','GCV','AIC','Risk')
  r_matrix <- matrix(0,nrow=4,ncol=T)
  rownames(r_matrix) <- c('CV','GCV','AIC','Risk')
  XY_matrix <- matrix(0,nrow=n,ncol=2*T)
  
  
  for(t in seq(1,T)) {
    set.seed(59087+j*t)
    # Normal error with mean 0 and var 1
    e <- rnorm(n, 0, 1)
    # sigma = 0.1
    e <- 0.1*e
    # Sample x 
    x <- sample_x(n,j)
    # Construct observerd y
    f <- true_f1(x)
    y <- f + e
    
    XY_matrix[,2*t-1] <- x
    XY_matrix[,2*t] <- y
    
    err_smallest <- unlist(optimize(function(z) { sum((f-H_lambda(x,z,p,nk)%*%y)^2) } ,c(0,100), tol=.Machine$double.eps))[2]
    
    opt_cv <- Optimal_Lambda(x,y,p,nk,'CV_score')
    opt_gcv <- Optimal_Lambda(x,y,p,nk,'GCV_score')
    opt_aicc <- Optimal_Lambda(x,y,p,nk,'AICc_score')
    opt_risk <- Optimal_Lambda(x,y,p,nk,'Risk_score')
    
    Lambda_matrix[1,t] <- opt_cv
    Lambda_matrix[2,t] <- opt_gcv
    Lambda_matrix[3,t] <- opt_aicc
    Lambda_matrix[4,t] <- opt_risk
    
    r_matrix[1,t] <- sum((f - H_lambda(x,opt_cv,p,nk)%*%y)^2)/err_smallest
    r_matrix[2,t] <- sum((f - H_lambda(x,opt_gcv,p,nk)%*%y)^2)/err_smallest
    r_matrix[3,t] <- sum((f - H_lambda(x,opt_aicc,p,nk)%*%y)^2)/err_smallest
    r_matrix[4,t] <- sum((f - H_lambda(x,opt_risk,p,nk)%*%y)^2)/err_smallest
    
  }
  
  r_matrix <- log(r_matrix)
  
  wilcox_rank <- vector(mode='list', length=4)
  names(wilcox_rank) <- c('CV', 'GCV', 'AICc', 'Risk')
  
  med_vec <- apply(r_matrix,1,median)
  
  rank_vec <- order(med_vec)
  
  test <- c()
  for(k in seq(1,3,1)) {
    r1 <- r_matrix[rank_vec[k],]
    r2 <- r_matrix[rank_vec[k+1],]
    ind <- (1-((r1-r2)==0)*1)==1
    if(sum(ind*1)==0) {
      test <- c(test,0)
    } else{
      test <- c(test,(wilcox.test(r1[ind],r2[ind],paired=T,alternative='less')$p.value < alpha)*1)
    }
  }
  
  if (sum(test) == 0) {
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- 2.5
    }
  } else if (sum(test) == 1) {
    pos <- which(test == 1)
    temp1 <- (1+pos)/2
    temp2 <- (pos+5)/2
    for(k in seq(1,pos)) {
      wilcox_rank[[rank_vec[k]]] <- temp1
    }
    for(k in seq(pos+1,4)) {
      wilcox_rank[[rank_vec[k]]] <- temp2
    }
  } else if (sum(test) == 2) {
    pos <- which(test == 0)
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- k
    }
    wilcox_rank[[rank_vec[pos]]] <- pos+0.5
    wilcox_rank[[rank_vec[pos+1]]] <- pos+0.5
  } else {
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- k
    }
  }
  
  print(unlist(wilcox_rank))
  
  # Plot true function
  plot(seq(0,1,1/100),true_f1(seq(0,1,1/100)),type='l',col='blue',main=paste('j=',j),xlab='',ylab='',ylim=c(-0.5,1.5))
  # Plot a typical simulation x,y
  points(XY_matrix[,1], XY_matrix[,2],col='red')
  
  boxplot_names <- c(wilcox_rank$CV,wilcox_rank$GCV,wilcox_rank$AICc,wilcox_rank$Risk)
  boxplot(r_matrix[1,],r_matrix[2,],r_matrix[3,],r_matrix[4,],names=boxplot_names,xlab='',ylab='',ylim=c(0,1.5))
  
  save(Lambda_matrix,XY_matrix,r_matrix,file=paste('design-density-J',j,'.RData',sep=''))
  
}
dev.off()
# ==========================Spatial Variation========================
jpeg('spatial-variation.jpg', width=4.25, height=3.25, units="in", res=1000, pointsize=4)
par(mfrow=c(3,4))

for(j in seq(1,J)) {
  
  Lambda_matrix <- matrix(0,nrow=4,ncol=T)
  rownames(Lambda_matrix) <- c('CV','GCV','AIC','Risk')
  r_matrix <- matrix(0,nrow=4,ncol=T)
  rownames(r_matrix) <- c('CV','GCV','AIC','Risk')
  XY_matrix <- matrix(0,nrow=n,ncol=2*T)
  
  
  for(t in seq(1,T)) {
    set.seed(12076+j*t)
    # Normal error with mean 0 and var 0.01
    e <- rnorm(n, 0, 1)
    # sigma = 0.1
    e <- 0.2*e
    # Sample x from uniform [0,1]
    x <- seq((1-0.5)/n,(n-0.5)/n,1/n)
    # Construct observerd y
    f <- true_f2(x,j)
    y <- f + e
    
    XY_matrix[,2*t-1] <- x
    XY_matrix[,2*t] <- y
    
    err_smallest <- unlist(optimize(function(z) { sum((f-H_lambda(x,z,p,nk)%*%y)^2) } ,c(0,100), tol=.Machine$double.eps))[2]
    
    opt_cv <- Optimal_Lambda(x,y,p,nk,'CV_score')
    opt_gcv <- Optimal_Lambda(x,y,p,nk,'GCV_score')
    opt_aicc <- Optimal_Lambda(x,y,p,nk,'AICc_score')
    opt_risk <- Optimal_Lambda(x,y,p,nk,'Risk_score')
    
    Lambda_matrix[1,t] <- opt_cv
    Lambda_matrix[2,t] <- opt_gcv
    Lambda_matrix[3,t] <- opt_aicc
    Lambda_matrix[4,t] <- opt_risk
    
    r_matrix[1,t] <- sum((f - H_lambda(x,opt_cv,p,nk)%*%y)^2)/err_smallest
    r_matrix[2,t] <- sum((f - H_lambda(x,opt_gcv,p,nk)%*%y)^2)/err_smallest
    r_matrix[3,t] <- sum((f - H_lambda(x,opt_aicc,p,nk)%*%y)^2)/err_smallest
    r_matrix[4,t] <- sum((f - H_lambda(x,opt_risk,p,nk)%*%y)^2)/err_smallest
    
  }
  
  r_matrix <- log(r_matrix)
  
  wilcox_rank <- vector(mode='list', length=4)
  names(wilcox_rank) <- c('CV', 'GCV', 'AICc', 'Risk')
  
  med_vec <- apply(r_matrix,1,median)
  
  rank_vec <- order(med_vec)
  
  test <- c()
  for(k in seq(1,3,1)) {
    r1 <- r_matrix[rank_vec[k],]
    r2 <- r_matrix[rank_vec[k+1],]
    ind <- (1-((r1-r2)==0)*1)==1
    if(sum(ind*1)==0) {
      test <- c(test,0)
    } else{
      test <- c(test,(wilcox.test(r1[ind],r2[ind],paired=T,alternative='less')$p.value < alpha)*1)
    }
  }
  
  if (sum(test) == 0) {
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- 2.5
    }
  } else if (sum(test) == 1) {
    pos <- which(test == 1)
    temp1 <- (1+pos)/2
    temp2 <- (pos+5)/2
    for(k in seq(1,pos)) {
      wilcox_rank[[rank_vec[k]]] <- temp1
    }
    for(k in seq(pos+1,4)) {
      wilcox_rank[[rank_vec[k]]] <- temp2
    }
  } else if (sum(test) == 2) {
    pos <- which(test == 0)
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- k
    }
    wilcox_rank[[rank_vec[pos]]] <- pos+0.5
    wilcox_rank[[rank_vec[pos+1]]] <- pos+0.5
  } else {
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- k
    }
  }
  
  print(unlist(wilcox_rank))
  
  # Plot true function
  plot(x,f,type='l',col='blue',main=paste('j=',j),xlab='',ylab='',ylim=c(-0.5,1.5))
  # Plot a typical simulation x,y
  points(XY_matrix[,1], XY_matrix[,2],col='red')
  
  boxplot_names <- c(wilcox_rank$CV,wilcox_rank$GCV,wilcox_rank$AICc,wilcox_rank$Risk)
  boxplot(r_matrix[1,],r_matrix[2,],r_matrix[3,],r_matrix[4,],names=boxplot_names,xlab='',ylab='',ylim=c(0,1.5))
  
  
  save(Lambda_matrix,XY_matrix,r_matrix,file=paste('spatial-variation-J',j,'.RData',sep=''))
  
}
dev.off()
# ==========================Variance Function========================
jpeg('variance-function.jpg', width=4.25, height=3.25, units="in", res=1000, pointsize=4)
par(mfrow=c(3,4))


for(j in seq(1,J)) {
  
  Lambda_matrix <- matrix(0,nrow=4,ncol=T)
  rownames(Lambda_matrix) <- c('CV','GCV','AIC','Risk')
  r_matrix <- matrix(0,nrow=4,ncol=T)
  rownames(r_matrix) <- c('CV','GCV','AIC','Risk')
  XY_matrix <- matrix(0,nrow=n,ncol=2*T)
  
  
  for(t in seq(1,T)) {
    set.seed(45097+j*t)
    
    # Sample x 
    x <- seq((1-0.5)/n,(n-0.5)/n,1/n)
    # Normal error with mean 0 and var 1
    e <- rnorm(n, 0, 1)
    # sigma 
    e <- 0.15*abs(1+0.4*(2*j-7)*(x-0.5))*e
    # Construct observerd y
    f <- true_f1(x)
    y <- f + e
    
    XY_matrix[,2*t-1] <- x
    XY_matrix[,2*t] <- y
    
    err_smallest <- unlist(optimize(function(z) { sum((f-H_lambda(x,z,p,nk)%*%y)^2) } ,c(0,100), tol=.Machine$double.eps))[2]
    
    opt_cv <- Optimal_Lambda(x,y,p,nk,'CV_score')
    opt_gcv <- Optimal_Lambda(x,y,p,nk,'GCV_score')
    opt_aicc <- Optimal_Lambda(x,y,p,nk,'AICc_score')
    opt_risk <- Optimal_Lambda(x,y,p,nk,'Risk_score')
    
    Lambda_matrix[1,t] <- opt_cv
    Lambda_matrix[2,t] <- opt_gcv
    Lambda_matrix[3,t] <- opt_aicc
    Lambda_matrix[4,t] <- opt_risk
    
    r_matrix[1,t] <- sum((f - H_lambda(x,opt_cv,p,nk)%*%y)^2)/err_smallest
    r_matrix[2,t] <- sum((f - H_lambda(x,opt_gcv,p,nk)%*%y)^2)/err_smallest
    r_matrix[3,t] <- sum((f - H_lambda(x,opt_aicc,p,nk)%*%y)^2)/err_smallest
    r_matrix[4,t] <- sum((f - H_lambda(x,opt_risk,p,nk)%*%y)^2)/err_smallest
    
  }
  
  r_matrix <- log(r_matrix)
  
  wilcox_rank <- vector(mode='list', length=4)
  names(wilcox_rank) <- c('CV', 'GCV', 'AICc', 'Risk')
  
  med_vec <- apply(r_matrix,1,median)
  
  rank_vec <- order(med_vec)
  
  test <- c()
  for(k in seq(1,3,1)) {
    r1 <- r_matrix[rank_vec[k],]
    r2 <- r_matrix[rank_vec[k+1],]
    ind <- (1-((r1-r2)==0)*1)==1
    if(sum(ind*1)==0) {
      test <- c(test,0)
    } else{
      test <- c(test,(wilcox.test(r1[ind],r2[ind],paired=T,alternative='less')$p.value < alpha)*1)
    }
  }
  
  if (sum(test) == 0) {
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- 2.5
    }
  } else if (sum(test) == 1) {
    pos <- which(test == 1)
    temp1 <- (1+pos)/2
    temp2 <- (pos+5)/2
    for(k in seq(1,pos)) {
      wilcox_rank[[rank_vec[k]]] <- temp1
    }
    for(k in seq(pos+1,4)) {
      wilcox_rank[[rank_vec[k]]] <- temp2
    }
  } else if (sum(test) == 2) {
    pos <- which(test == 0)
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- k
    }
    wilcox_rank[[rank_vec[pos]]] <- pos+0.5
    wilcox_rank[[rank_vec[pos+1]]] <- pos+0.5
  } else {
    for(k in seq(1,4,1)) {
      wilcox_rank[[rank_vec[k]]] <- k
    }
  }
  
  print(unlist(wilcox_rank))
  
  # Plot true function
  plot(x,f,type='l',col='blue',main=paste('j=',j),xlab='',ylab='',ylim=c(-0.5,1.5))
  # Plot a typical simulation x,y
  points(XY_matrix[,1], XY_matrix[,2],col='red')
  
  boxplot_names <- c(wilcox_rank$CV,wilcox_rank$GCV,wilcox_rank$AICc,wilcox_rank$Risk)
  boxplot(r_matrix[1,],r_matrix[2,],r_matrix[3,],r_matrix[4,],names=boxplot_names,xlab='',ylab='',ylim=c(0,1.5))
  
  
  save(Lambda_matrix,XY_matrix,r_matrix,file=paste('variance-function-J',j,'.RData',sep=''))
  
}
dev.off()


