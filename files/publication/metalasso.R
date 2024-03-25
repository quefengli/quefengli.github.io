library(penalized)
library(glmnet)

## calculate the optimization problem with a single lam
## Args:
##      X.all:     all X covariates
##      Y.all:     all Y responses
##      obs:       number of observations in each dataset
##      lam:       penality of theta's
##      method:    either "penalized" or "glmnet"
##                 ("glemnet" could hand unbalanced samples size)
##      maxit:     maximal number of iterations
##      tol:       tolerance level
## Returns:
##      coe:       fitted coefficients
##      gamma:     fitted gamma's
##      iteration: number of iterations
##      converge:  indicator if convergence is achieved
##      lam:       penality used 
##      diff:      last step difference

metalasso <- function(X.all, Y.all, obs, lam1, method, maxit, tol){
  ## starting and ending index of each dataset
  start.idx   <- cumsum(obs) + 1 - obs    # starting index of each dataset
  end.idx     <- cumsum(obs)              # ending index of each dataset
  M           <- length(start.idx)        # number of datasets
  p           <- ncol(X.all)              # number of covariates
  N           <- sum(obs)                 # total number of obserations
  gamma       <- matrix(NA, p, maxit + 1) # iterations of gamma
  gamma[, 1]  <- rep(1, p)
  theta       <- matrix(NA, M * p, maxit) # iterations of theta
  X.tha       <- matrix(NA, N, p)         # colMultiply(X.all, theta)
  beta.hat    <- vector("list", M)        # iterations of beta.hat
  coef        <- vector("list", M)        # final estimate of coefficients
  itr         <- 1
  m.diff      <- NULL                     # marginal error
  for (m in 1:M){
    beta.hat[[m]] <- matrix(NA, p, maxit)
  }
    
  if (method == "glmnet") {
    {
    while(!(itr > maxit)){
      ## Iterate as: theta --> gamma --> theta
      
      for (m in 1:M){
      ## In each dataset, fit Y.all ~ colMultiply(X.all, gamma*rho)
        theta.fit <- glmnet(t(t(X.all[start.idx[m]:end.idx[m], ]) * 
                              gamma[, itr]),
                            Y.all[start.idx[m]:end.idx[m]],
                            alpha = 1, 
                            family = "binomial", 
                            lambda = lam1,
                            standardize = FALSE)
        
        theta[((m - 1) * p + 1):(m * p), itr] <- as.vector(theta.fit$beta)
        
        ## adjust X.all by colMultiply(X.all, theta) for further usage
        X.tha[start.idx[m]:end.idx[m], ] <- t(t(X.all[start.idx[m]:end.idx[m], ]) *
                                            as.vector(theta.fit$beta))

        beta.hat[[m]][, itr] <- theta[((m - 1) * p + 1):(m * p), itr] * gamma[, itr]
        ## calculate iteration difference
        if (itr == 1) {
          m.diff[m] <- max(abs(beta.hat[[m]][, itr]))
        }
        else {
          m.diff[m] <- max(abs(beta.hat[[m]][, itr] - beta.hat[[m]][, itr - 1]))
        }
      }
      
      if(max(m.diff) < tol) break         # break iterations if diff < tol

      itr <- itr + 1
            
      gamma.fit <- glmnet(X.tha, Y.all, alpha = 1, family = "binomial", 
                          lambda = 1, 
                          weights = rep(1 / obs, obs), 
                          standardize = FALSE)
      gamma[, itr] <- as.vector(gamma.fit$beta)
    }
  }
  }
    
  if (method == "penalized") {
    while(!(itr > maxit)){
      ## Iterate as: theta --> gamma --> theta
            
      for (m in 1:M){
      ## In each dataset, fit Y.all ~ colMultiply(X.all, gamma*rho)
        theta.fit <- penalized(Y.all[start.idx[m]:end.idx[m]],
                               t(t(X.all[start.idx[m]:end.idx[m], ]) *
                                 gamma[, itr]),
                               positive = TRUE, 
                               unpenalized = ~0, model = "logistic",
                               lambda1 = lam1,
                               standardize = FALSE, trace = FALSE)
        theta[((m - 1) * p + 1):(m * p), itr] <- theta.fit@penalized
        
        ## adjust X.all by colMultiply(X.all, theta) for further usage
        X.tha[start.idx[m]:end.idx[m], ] <- t(t(X.all[start.idx[m]:end.idx[m], ]) *
                                            theta.fit@penalized)
        
        beta.hat[[m]][, itr] <- theta[((m - 1) * p + 1):(m * p), itr] * gamma[, itr]
        ## calculate iteration difference
        if (itr == 1) {
          m.diff[m] <- max(abs(beta.hat[[m]][, itr]))
        }
        else {
          m.diff[m] <- max(abs(beta.hat[[m]][, itr] - beta.hat[[m]][, itr - 1]))
        }
      }
      
      if(max(m.diff) < tol) break         # break iterations if diff < tol
      
      itr <- itr + 1

      gamma.fit <- penalized(Y.all, X.tha, 
                             unpenalized = ~0, model = "logistic",
                             lambda1 = 0.005,
                             standardize = FALSE, trace = FALSE)
      gamma[, itr] <- gamma.fit@penalized
    }
  }


  ## determine if convergence is achieved
  if (itr == 1) {
    iteration <- itr
    converge  <- FALSE
  }
  else {
    if (itr > maxit) {
      iteration <- itr - 1
      converge  <- FALSE
    }
    else {
      iteration <- itr
      converge  <- TRUE
    }  
  }

  for (m in 1:M){
    coef[[m]] <- beta.hat[[m]][, iteration]
  }
  
  return(list(coe         = coef,
              gamma       = gamma[, iteration], 
              iteration   = iteration,
              converge    = converge,
              lam         = lam1,
              diff        = max(m.diff)
              ))
}
