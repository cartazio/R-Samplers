# Batch VB for LDA
library(MCMCpack)

# Inputs: n_iv, K, alpha_k, gamma_v

dat <- matrix(rpois(50, lambda = 1), nrow = 5)
colnames(dat) <- letters[1:10]
dat

# n_iv = number of times word v occurs in doc i (observed)
# DOCUMENT TERM MATRIX
# i rows, v columns: n[i,v]
n <- dat

# K number of topics, set my researcher
K <- 3

# Hyperparameters set by researcher
# one alpha for each K
# hyperparameter for pi (all n pi's have same alpha's)
alpha <- rep(1, times = K)

# example: individual distributions over topics
rdirichlet(nrow(n), alpha)

# one gamma for each v
# hyperparameter for B (all b's have same gamma's)
gamma <- rep(1, times = ncol(n))

# example: topic distributions over words
rdirichlet(K, gamma)

# set b_vk DELETE
#b_vp <- matrix(rep(gamma, times = K), nrow = ncol(n), ncol = K)

# array(rep(n, times = 3), dim = c(5, 10, 3))

# expectation step
estep <- function(n, b_vp, alpha) {
    # update pi's and c's
    # Initialize pi_ik_vp: pi[i,k] variational parameter
    pi_vp <- rdirichlet(nrow(n), alpha)
    # Empty c_ivk_vp: c[i,v,k] variational parameter array
    c_vp <- array(0, dim = c(nrow(n), ncol(n), K))
    
    #iter <- 0
    #repeat {
      
      #iter <- iter + 1
      
      for (i in 1:nrow(n)) {
        
      pi_vp_old <- pi_vp
      pi_vp <- matrix(rep(alpha, times = nrow(n)), nrow = nrow(n), ncol = K)
        
      # for v words
      for (v in 1:ncol(n)) {
        
        # for k topics
        for(k in 1:K) {
          # Added column indexes here, is that ok?
          c_vp[i,v,k] <- exp(digamma(b_vp[v,k]) + digamma(pi_vp_old[i,k])) #digammas per K
          }
        c_vp[i,v,] <- c_vp[i,v,]/sum(c_vp[i,v,]) #this is normalize to a simplex?
        pi_vp[i,k] <- pi_vp[i,k] + (n[i,v] * c_vp[i,v,k]) #what kind of multiplication is this? 
        }
      #if (iter > 10){ # (1/K)*sum(abs(pi_vp[i,k] - pi_vp_old[i,k])) < thresh) {
        #break
    }
      #}
    #}
    list(pi_vp, c_vp)
}

#
# estep(n, matrix(rep(gamma, times = K), nrow = ncol(n), ncol = K), alpha)
#

# # maximization step
# mstep <- function (gamma, s) {
#   # b_vp <- matrix(0, nrow = ncol(n), ncol = K) # this is initialized in lda_vi so dont need it here
#   for (v in 1:ncol(n)) { #added in some extra indexes, is this right?
#   for (k in 1:K) {
#     b_vp[v,k] <- gamma[v] + s[v,k]
#   }
#   }
#   b_vp
# }
# 
# #
# mstep(gamma, matrix(0, nrow = ncol(n), ncol = K))
# #

# optimize
lda_vi <- function (n, K, alpha, gamma) {
  # estimate b_vk_vp using EM for multinomial mixtures
  b_vp <- t(rdirichlet(K, gamma)) # fill this in with random init
  
  # initialize counts n_iv
  n <- dat # document term matrix with counts

  iter = 0
  while (iter < 10) {
    
    iter = iter + 1
    s = matrix(0, nrow = ncol(n), ncol = K) # expected sufficient statistic
    
    # e step
    #for (i in 1:nrow(n)) {
      estep_out <- estep(n, b_vp, alpha)
      pi_vp <- estep_out[[1]]
      c_vp <- estep_out[[2]]
      for (i in 1:nrow(n)) {
      for (k in 1:K) {
        for (v in i:ncol(n)) {
          s[v,k] = s[v,k] + (n[i,v] * c_vp[i,v,k]) # what kind of multiplication is this?
      }
      }
    }
    pi_vp
    c_vp
    s
    # m step
    for (v in 1:ncol(n)) { #added in some extra indexes, is this right?
        for (k in 1:K) {
          b_vp[v,k] <- gamma[v] + s[v,k]
        }
    }
    b_vp
    #if (iter < 10) {converged <- TRUE} # need to compute elbo and log likelihood
  }
  list(pi_vp, b_vp, c_vp, s, K)
  }



lda_vi(n, K, alpha, gamma)  
  
  
  
  
  
  
  
  


