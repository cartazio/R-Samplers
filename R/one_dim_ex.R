## Gibbs sampler

# x an n vector of data
# pi a k vector
# mu a k vector
sample_z = function(x,pi,mu){
  dist_matrix = outer(mu, x, "-")
  # distance matrix, x from mu
  # k by n matrix, d_kj =(mu_k - x_j)
  p.z = log(as.vector(pi)) + dnorm(dist_matrix,0,1, log = TRUE) # sample z given pi and distance to mu ON LOG SCALE
  p.z = apply(p.z, 2, function(x){return(x/sum(x))}) # normalize columns
  
  z = rep(0, length(x))
  for(i in 1:length(z)){
    z[i] = sample(1:length(pi), size = 1, prob = p.z[,i], replace = TRUE) #Should be cat/multinomial?
  }
  return(z)
}

#sample_z(x, pi, mu)

# z an n vector of cluster allocations (1...k)
# k the number of clusters
sample_pi = function(z,k){
  counts = colSums(outer(z,1:k,FUN="=="))
  pi = MCMCpack::rdirichlet(1,counts+1)
  return(pi)
}

#sample_pi(z, k)

# x an n vector of data
# z an n vector of cluster allocations
# k the number o clusters
# prior.mean the prior mean for mu
# prior.lambda the prior precision for mu
sample_mu = function(x, z, k, prior){
  df = data.frame(x = x, z = z)
  mu = rep(0, k)
  for(i in 1:k){
    n_k = sum(z == i) #number of data assigned to component k
    x_bar = ifelse(n_k == 0, 0, mean(x[z == i])) #kth component's mean
    
    lambda_hat = 1/(n_k + prior$lambda) #lambda is the prior sd on the mixture components
    
    mu_hat = log(x_bar) + (log(n_k) - log(n_k + 1))
      #log(prior$mean) + log(prior$lambda) + log(x_bar) + log(n_k) - log(lambda_hat)
    
    mu[i] = rnorm(1, mean = mu_hat, sd = sqrt(lambda_hat))
  }
  return(mu)
}

#sample_mu(x, z, k, prior)

Gauss1dim <- function(x, k){
  gibbsHarness(
    ,InitGibbState=function(k, x){
      mu = rnorm(k,0,1)
      z = sample_z(x,pi,mu)
      }
    
    ,TransitionProposal=function(x){rnorm(1) }# doesn't care about current state
    
    
    ,ApplyTransition=function(state,proposal){state + proposal }
    
    ,ShouldWeTerminate=function(step,state,proposal){ (step > 10 )
    }
  )
  
}
