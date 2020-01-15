## Gibbs sampler
library(tidyverse)

source("./R/GibbsCore.R")

# x is an n vector of data
# pie is a k vector
# mu is a k vector
sample_z = function(x, pie, mu){
  k = length(pie)
  dist_matrix = outer(mu, x, "-")
  
  # distance matrix, x from mu
  # k by n matrix, d_kj = (mu_k - x_j)
  prob_z = exp(log(as.vector(pie)) + log(dnorm(dist_matrix, 0, 1))) # sample z given pie and distance to mu ON LOG SCALE
  prob_z = apply(prob_z, 2, function(x){return(x/sum(x))}) # normalize columns
  
  #z = rep(0, length(x))
  z = matrix(0, nrow = nrow(prob_z), ncol = ncol(prob_z))

    for(i in 1:ncol(z)){
    #z[i] = sample(1:length(pie), size = 1, replace = TRUE, prob = prob_z[,i]) #Should be cat/multinomial?
    z[,i] <- rmultinom(1, 1, prob = prob_z[,i])
  }
  return(z)
}

# z is a matrix of cluster allocations (n*k)
# k is the number of clusters
sample_pie = function(z, k){
  #counts = colSums(z)
  pie = MCMCpack::rdirichlet(1,rep(1,k))
  return(pie)
}

# x is an n vector of data
# z is a matrix of cluster allocations
# k is the number o clusters
# lambda is the prior sd for mu
sample_mu = function(x, z, k, prior=list(mean=0,prec=0.1)){
  mu = rep(0, k)
  for(i in 1:k){
    n_k = colSums(z) # number of data assigned to component k
    x_bar = colSums(x * z) # kth component's mean
    
    # lambda_hat = log(1) - log(n_k + lambda) #lambda is the prior sd on the mixture components
    # mu_hat = log(x_bar) + (log(n_k) - log(n_k + 1))
    
    lambda_hat = n_k + prior$prec
    mu_hat = (prior$mean * prior$prec + x_bar * n_k)/lambda_hat ## Not logged but doesnt work when logged...
    
    mu[i] = rnorm(1, mean = mu_hat, sd = sqrt(1/lambda_hat))
  }
  return(mu)
}

Gauss1dim <- function(k, x){
  gibbsLoop(
    InitGibbState=function(){
      mu = rnorm(k, 0, 1)
      pie = rep(1/k, k)
      z = sample_z(x, pie, mu)
      list(mu = mu, pie = pie, z = z)
      }
    
    ,TransitionProposal=function(previousState){ # p_* for proposed_*
     p_mu <- sample_mu(x, previousState$z, k)
     p_pie <- sample_pie(previousState$z, k)
     p_z <- sample_z(x, p_pie, p_mu)
     list(mu = p_mu, pie = p_pie, z = p_z)
    }
    
    ,ApplyTransition=function(previousState,proposal){
      mu <- proposal$mu
      pie <- proposal$pie
      z <- proposal$z
      list(mu = mu, pie = pie, z = z)
      }
    
    ,ShouldWeTerminate=function(step,state,proposal){(step > 100)
    }
  )
}

