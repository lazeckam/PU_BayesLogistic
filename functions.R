

sigma <- function(t){
  
  return(1/(1 + exp(-t)))
  
}


stepPUPG <- function(S, X,
                     beta_b, beta_B, beta_start, 
                     c_alpha, c_beta, c_start,
                     beta_B_inv=NULL, beta_B_inv_beta_b=NULL){
  # one step of Gibbs sampler for PU data
  
  if(is.null(beta_B_inv)){
    beta_B_inv <- solve(beta_B)
  }
  if(is.null(beta_B_inv_beta_b)){
    beta_B_inv_beta_b <- beta_B_inv %*% beta_b
  }
  
  # sampling p(omega|beta)
  omega_new <- sapply(X %*% beta_start, BayesLogit::rpg.devroye, num=1, h=1)
  
  # sampling p(y|beta, s, c)
  #  - sampling p(y|beta)
  sigma_new <- sigma(X %*% beta_start)
  #  - multiplying p(y|beta) by p(s|y, c)
  pY1 <- sigma_new * ifelse(S == 1, c_start, 1 - c_start)
  pY0 <- (1 - sigma_new) * ifelse(S == 1, 0, 1)
  pY1_new <- pY1/(pY1 + pY0)
  #  - sampling y
  Y_new <- sapply(pY1_new, rbinom, n=1, size=1)
  
  # sampling p(beta|omega, y)
  #  - computing posterior parameters
  beta_B_post <- chol2inv(chol(t(X) %*% diag(omega_new) %*% X + beta_B_inv))
  beta_b_post <- beta_B_post %*% (t(X) %*% (Y_new - 1/2) + beta_B_inv_beta_b)
  #  - sampling beta
  beta_new <- as.vector(mvtnorm::rmvnorm(1, beta_b_post, beta_B_post))
  
  # sampling p(c|s, y)
  S_factor <- factor(S, levels=0:1)
  Y_new_factor <- factor(Y_new, levels=0:1)
  #  - computing posterior parameters
  tabSY_new <- table(S_factor, Y_new_factor)
  c_alpha_post <- c_alpha + tabSY_new[2,2]
  c_beta_post <- c_beta + tabSY_new[1,2]
  #  - sampling c
  c_new <- rbeta(1, c_alpha_post, c_beta_post)
  
  results <- list(beta =beta_new,
                  c    =c_new,
                  omega=omega_new,
                  Y    =Y_new,
                  pY1  =pY1_new)
  
  return(results)
}


PUPG <- function(S, X, # data - an n-dimensional vector and nxp-dimensional matrix
                 N, B, # N - number of required samples, B - number of samples to discard from the beginning of the chain ("burn-out")
                 beta_b=rep(0, ncol(X)), beta_B=diag(rep(100, ncol(X))), beta_start=beta_b, # hyperparameters & initial value for beta
                 c_alpha=1, c_beta=1, c_start=0.5){ # hyperparameters & initial value for C
  
  p <- ncol(X)
  n <- nrow(X)
  
  beta_B_inv <- chol2inv(chol(beta_B))
  beta_B_inv_beta_b <- beta_B_inv %*% beta_b
  
  results <- list(beta_values =array(NA, dim=c(N, p), dimnames=list(N=1:N, coefficients=paste0("coef", 0:(p-1)))),
                  c_values    =array(NA, dim=N, dimnames=list(N=1:N)),
                  omega_values=array(NA, dim=c(N, n), dimnames=list(N=1:N, n=1:n)),
                  Y_values    =array(NA, dim=c(N, n), dimnames=list(N=1:N, n=1:n)),
                  pY1_values  =array(NA, dim=c(N, n), dimnames=list(N=1:N, n=1:n)))
  
  results[["info"]] <- list(N=N, B=B,
                            beta_b=beta_b, beta_B=beta_B, beta_start=beta_start,
                            c_alpha=c_alpha, c_beta=c_beta, c_start=c_start)
  
  beta_old <- beta_start
  c_old <- c_start
  
  for(i in 1:(B + N)){
    
    res_tmp <- stepPUPG(S=S, X=X,
                        beta_b=beta_b, beta_B=beta_B, beta_start=beta_old, 
                        c_alpha=c_alpha, c_beta=c_beta, c_start=c_old,
                        beta_B_inv=beta_B_inv, beta_B_inv_beta_b=beta_B_inv_beta_b)
    
    beta_old <- res_tmp$beta
    c_old <- res_tmp$c
    
    if(i > B){
      results[["beta_values"]][i - B,] <- res_tmp$beta
      results[["c_values"]][i - B] <- res_tmp$c
      results[["omega_values"]][i - B,] <- res_tmp$omega
      results[["Y_values"]][i - B,] <- res_tmp$Y
      results[["pY1_values"]][i - B,] <- res_tmp$pY1
    }
  }
  
  return(results)
}