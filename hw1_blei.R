library(tidyverse)
library(janitor)
#install.packages("MCMCpack")
library(MCMCpack)

source("./R/GibbsCore.R")

## Read in data
undata <- read_csv("./data/country_profiles.csv") %>% 
  clean_names() %>% 
  mutate_at(vars(3:50), as.numeric) %>% 
  na_if(-99)

# Initial state for mu's is random
# Sample z's from x's, mu's, theta (set)
# Sample mu's from x's and z's
# Repeat until converged

fakedat <- matrix(1:150, ncol = 10)

mysampler <- function(clustercount=4
                      ,mydata=fakedat
                      ,dirichlet_alpha=c(1/4,1/4,1/4,1/4)
                      ,LAMBDA=1) {
  
   # need to check sum to 1. abort otherwise 
  clustercount <- 4
  data_dim <- ncol(fakedat)
  rows <- nrow(fakedat)
  mu_min <- apply(fakedat, 2, min)
  mu_max <- apply(fakedat, 2, max)
  mu_median <- apply(fakedat, 2, median)
  mu_sd <- apply(fakedat, 2, sd)
  

  logSimpleMVNPDF <- function(data,mu) {
    logProd=0 # log 1 == 0
    for (i in 1:length(mu)) {
      logProd = logProd + dnorm(data[i], mean = mu[i], sd = LAMBDA, log = TRUE)
    }

    logProd
  }

  logVectored_MVN_PDF <- function(data,muMultiVector){
    #  each row is a mu for a cluster
    # for each  column J of row K , we calculation prob x_J is from Normal(U_{J,K} , LAMBDA)
    # so the mvn pdf for data at a row would be  Product_j  {prob(x_{j} in pdfNorm(u_{j,k},lambda)}
    result <- apply(muMultiVector,1,function(mu) {logSimpleMVNPDF(data,mu)})
    print(result)
    result
  }
  ## the >= 0.001 clamping should be a hyper parameter
  ## this is to ensure every cluster is non empty with nonzero probability
  clamped_theta_gen <- function(){ map2_dbl(  rdirichlet(1, dirichlet_alpha), rep(0.001, clustercount),function(x,y){max(x,y)})}
  theta_init <-clamped_theta_gen()
  
  
  # Initiate mu's, draw from known mean and sd
  init_mu <- function(dummy){rnorm(data_dim, mean = mu_median, sd = mu_sd)}
  zero_matrix <- matrix(0, nrow = clustercount, ncol = data_dim)
  
  
  # for a cluster_count hyper param, assignment that isn't nonempty for every cluster is invalid, try again!
  z_rejection_sampler_init <- function(theta){
     z_guess <- matrix(0, nrow = rows, ncol = clustercount)
      while(!(  every(apply(z_guess,2,sum),function(x){x > 0})) 
        # ||  !(sum(z_guess)-rows) < 0.001 
        )  {z_guess <- t(rmultinom(rows,1,theta))
  }
     z_guess
  }
  

  z_rejection_sampler_clusterWeights <-function(memberProbWeights){
    print(memberProbWeights)
     z_guess <- matrix(0, nrow = rows, ncol = clustercount)
      while(!(  every(apply(z_guess,2,sum),function(x){x > 0})) ){

        z_guess <- t(apply(memberProbWeights,1,function(row_Theta){
           print(row_Theta)
            t(rmultinom(1,1,row_Theta))

        } ))
      } 

  }
  
  gibbsHarness(
    InitGibbState= function(){
      list(
         # THETA = theta_init# MU= t(apply(zero_matrix, 1, init_mu)  )  # initialize all mu's, k x
        Z = z_rejection_sampler_init(theta_init)  # initialize the cluster assignment Z
      )}
    
    ,TransitionProposal=
      function(previousState){
        # pMu <- previousState$MU
        pZ <- previousState$Z 
        newTheta <- rdirichlet(1, dirichlet_alpha) 
        cluster_sum <- t(matrix(0, nrow = clustercount, ncol = data_dim))
          ZeeSUM <- apply(pZ, 2, sum)
          for (n1 in 1:clustercount) { 
            for (r in 1:rows ) {
              if (pZ[r,n1]== 1) {
                 cluster_sum[,n1]= fakedat[r,] + cluster_sum[,n1] #fix in the morning
                } } } 
        cluster_ave <- cluster_sum/ZeeSUM
              
        list(THETA= clamped_theta_gen() , MU= exp(log(cluster_ave) + log(ZeeSUM) - log(1+ ZeeSUM) ) )
              
      } 
    
    ##
    ,ApplyTransition=function(ignorestate,proposal){
        proTHETA <- proposal$THETA 
        proMU <- proposal$MU 
        log_prob_matrix <- apply(mydata, 1, function(myrow){logVectored_MVN_PDF(myrow,proMU)})
        
          list(Z= z_rejection_sampler_clusterWeights(exp(  abs(max(log_prob_matrix)) + 10 + log(proTHETA)+ 
                  log_prob_matrix
                        
                    ))
             )
          


    }
    ,ShouldWeTerminate=function(step,state,proposal){(step > (rows *  clustercount * data_dim)^2 )}
      # estimate invented via "Tuchas" based estimation  
  )   
}

mysampler()
