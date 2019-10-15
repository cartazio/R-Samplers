## Gibbs sampler

# x is an n vector of data
# pi is a k vector
# mu is a k vector
sample_z = function(x, pi, mu){
  dist_matrix = outer(mu, x, "-")
  # distance matrix, x from mu
  # k by n matrix, d_kj = (mu_k - x_j)
  prob_z = log(as.vector(pi)) + log(dnorm(dist_matrix, 0, 1)) # sample z given pi and distance to mu ON LOG SCALE
  prob_z = apply(prob_z, 2, function(x){return(x/sum(x))}) # normalize columns
  
  #z = rep(0, length(x))
  z = matrix(0, nrow = length(x), ncol = k)
    
  for(i in 1:nrow(z)){
    #z[i] = sample(1:length(pi), size = 1, replace = TRUE, prob = prob_z[,i]) #Should be cat/multinomial?
    z[i,] <- rmultinom(1,1,prob = prob_z[,i])
  }
  
  return(z)
}

#sample_z(x, pi, mu)

# z is an n vector of cluster allocations (1...k)
# k is the number of clusters
sample_pi = function(z, k){
  counts = colSums(outer(z, 1:k, FUN="=="))
  pi = MCMCpack::rdirichlet(1,counts+1)
  return(pi)
}

#sample_pi(z, k)

# x is an n vector of data
# z is an n vector of cluster allocations
# k is the number o clusters
# lambda the prior sd for mu
sample_mu = function(x, z, k, prior=list(mean=0,prec=0.1)){
  df = data.frame(x = x, z = z)
  mu = rep(0, k)
  for(i in 1:k){
    n_k = sum(z == i) #number of data assigned to component k
    x_bar = ifelse(n_k == 0, 0, mean(x[z == i])) #kth component's mean
    
    # lambda_hat = log(1) - log(n_k + lambda) #lambda is the prior sd on the mixture components
    # mu_hat = log(x_bar) + (log(n_k) - log(n_k + 1))
    
    lambda_hat = n_k+prior$prec
    mu_hat = (prior$mean * prior$prec + x_bar * n_k)/lambda_hat ## Not logged but doesnt work when logged...
    
    mu[i] = rnorm(1, mean = mu_hat, sd = sqrt(1/lambda_hat))
  }
  return(mu)
}

#sample_mu(x, z, k)

Gauss1dim <- function(k, x){
  gibbsHarness(
    InitGibbState=function(){
      mu = rnorm(k, 0, 1)
      pi = rep(1/k, k)
      z = sample_z(x, pi, mu)
      list(mu = mu, pi = pi, z = z)
      }
    
    ,TransitionProposal=function(previousState){ # p_* for proposed_*
     p_mu <- sample_mu(x, previousState$z, k)
     p_pi <- sample_pi(previousState$z, k)
     p_z <- sample_z(x, p_pi, p_mu)
     list(mu = p_mu, pi = pi, z = p_z)
    }
    
    ,ApplyTransition=function(previousState,proposal){
      mu <- proposal$mu
      pi <- proposal$pi
      z <- proposal$z
      list(mu = mu, pi = pi, z = z)
      }
    
    ,ShouldWeTerminate=function(step,state,proposal){(step > 10)
    }
  )
}

x <- fakedat[,1]
Gauss1dim(k = 1, x = x)
