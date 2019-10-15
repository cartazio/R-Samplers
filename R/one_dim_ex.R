## Gibbs sampler

# x an n vector of data
# pi a k vector
# mu a k vector
sample_z = function(x, pi, mu){
  dist_matrix = outer(mu, x, "-")
  # distance matrix, x from mu
  # k by n matrix, d_kj =(mu_k - x_j)
  prob_z = log(as.vector(pi)) + dnorm(dist_matrix, 0, 1, log = TRUE) # sample z given pi and distance to mu ON LOG SCALE
  prob_z = apply(prob_z, 2, function(x){return(x/sum(x))}) # normalize columns
  
  z = rep(0, length(x))
  for(i in 1:length(z)){
    z[i] = sample(1:length(pi), size = 1, replace = TRUE, prob = prob_z[,i]) #Should be cat/multinomial?
  }
  return(z)
}

#sample_z(x, pi, p_mu)

# z an n vector of cluster allocations (1...k)
# k the number of clusters
sample_pi = function(z, k){
  counts = colSums(outer(z, 1:k, FUN="=="))
  pi = MCMCpack::rdirichlet(1,counts+1)
  return(pi)
}

#sample_pi(z, k)

# x an n vector of data
# z an n vector of cluster allocations
# k the number o clusters
# prior.mean the prior mean for mu
# prior.lambda the prior precision for mu
sample_mu = function(x, z, k, lambda = 1){
  df = data.frame(x = x, z = z)
  mu = rep(0, k)
  for(i in 1:k){
    n_k = sum(z == i) #number of data assigned to component k
    x_bar = ifelse(n_k == 0, 0, mean(x[z == i])) #kth component's mean
    
    lambda_hat = log(1) - log(n_k + lambda) #lambda is the prior sd on the mixture components
    
    mu_hat = log(x_bar) + (log(n_k) - log(n_k + 1))
      #log(prior$mean) + log(prior$lambda) + log(x_bar) + log(n_k) - log(lambda_hat)
    
    mu[i] = rnorm(1, mean = mu_hat, sd = lambda_hat)
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
      }
    
    ,TransitionProposal=function(mu, pi){ # p_* for proposed_*
     p_mu <- sample_mu(x, z, k)
     p_pi <- sample_pi(z, k)
     p_z <- sample_z(x, pi, p_mu)
    }
    
    ,ApplyTransition=function(state,proposal){
      mu <- p_mu
      pi <- p_pi
      z <- p_z
      }
    
    ,ShouldWeTerminate=function(step,state,proposal){(step > 10)
    }
  )
}

Gauss1dim(k = 1, x = x)
