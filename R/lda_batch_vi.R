# Batch VB for LDA
library(MCMCpack)
library(abind)

# Inputs: 
# n_iv = document-term matrix
# K = number of topics, set by researcher

# Hyperparameters set by researcher

# one alpha for each K
# hyperparameter for pi (all n pi's have same alpha's)

# one gamma for each v
# hyperparameter for B (all b's have same gamma's)

# expectation step
estep <- function(n, b_vp ,alpha,ix) {
     K <- length(alpha)
    # update pi's and c's
    # Initialize pi_ik_vp: pi[i][k] variational parameter, i fixed
    pi_vp <- alpha
    # Empty c_ivk_vp: c[i][v,k] variational parameter array, i fixed
    c_vp <- array(1/K, dim = c(ncol(n), K))
    
    iter <- 0
    elbo_i <- c()
    repeat {
      
      iter <- iter + 1
        
      pi_vp_old <- pi_vp
      pi_vp <- alpha
        
      # for v words
      for (v in 1:ncol(n)) {
        
        # for k topics
        for(k in 1:K) {
            c_vp[v,k] <- exp(digamma(b_vp[v,k]) + digamma(pi_vp_old[k])) # digammas per K
          }
        c_vp[v,] <- c_vp[v,]/sum(c_vp[v,]) #this is normalized to a simplex
        pi_vp <- pi_vp + (n[ix,v] * c_vp[v,]) 
      }
      elbo_i <- sum((pi_vp - pi_vp_old)^2)
      if (iter > 10){ 
        break
    }
      }
    list(pi_vp = pi_vp, c_vp = c_vp, pi_vp_old = pi_vp_old, elbo_i = elbo_i)
}

# optimize
lda_vi <- function (n, K, alpha, gamma, max_iter) {
  b_vp <- t(rdirichlet(K, gamma)) # fill this in with random init

  iter = 0
  pi_ik <- matrix(0, nrow=nrow(n), ncol = K)
  s = matrix(0, nrow = ncol(n), ncol = K) # expected sufficient statistic
  
  trace_pi <- c()
  trace_b <- c()
  elbo <- c()
  while (iter < max_iter) {
    
    iter = iter + 1
 
    # e step
    for (i in 1:nrow(n)) {
      estep_out <- estep(n, b_vp, alpha, i)
      pi_vp <- estep_out[[1]]
      c_vp <- estep_out[[2]]
      pi_vp_old <- estep_out[[3]]
      elbo_i <- estep_out[[4]]
      
      pi_ik[i,]<- pi_vp
      s = s + (n[i,] * c_vp[,]) 
      s
      }
    pi_ik
    
    # m step
    for (v in 1:ncol(n)) {
        for (k in 1:K) {
          b_vp[v,k] <- gamma[v] + s[v,k]
        }
    }
    # Calculate ELBO = E[log(joint)] - E[log(entropy)]
    # elbo_i <- elbo_i + sum((pi_vp  - pi_vp_old)^2)
    elbo <- c(elbo, elbo_i)
    
    trace_pi <- abind(trace_pi, pi_ik, along=3)
    trace_b <- abind(trace_b, b_vp, along = 3)
 
    }
    list(pi_ik = pi_ik, b_vp = b_vp, s = s, K = K, trace_pi = trace_pi, 
         trace_b = trace_b, elbo = elbo, iter = 1:max_iter) # dont save c, too much space
  }


  
  
  
  
  
  
  


