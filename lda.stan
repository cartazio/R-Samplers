data {
int <lower=1> K; // num topics
int <lower=1> V; // num words
int <lower=1> M; // num docs
int <lower=0> X[M,V]; // doc-term matrix

vector <lower=0> [K] alpha; // hyperparam for theta prior
vector <lower=0> [V] beta; // hyperparam for phi prior
}

parameters {
simplex[K] theta[M]; // topic dist for doc m, MxK matrix
simplex[V] phi[K]; // word dist for topic k, KxV matrix
}

transformed parameters  {
  
  
}

model {
for (m in 1:M)
theta[m] ~ dirichlet(alpha); // prior

for (k in 1:K)
phi[k] ~ dirichlet(beta); // prior

for (m in 1:M) {
  for (v in 1:V) {
    real gamma[K];
      for (k in 1:K) {
        gamma[k] = log(theta[m,k]) + log(phi[k,v]); // prob(word | theta, phi) = theta * phi 
      }
    target += log_sum_exp(gamma); // likelihood;
  }
}
}

