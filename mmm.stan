data {
int<lower=1> K; // num topics
int<lower=1> V; // num words
int<lower=1> M; // num docs
matrix <lower=0> [M,V] X; // doc-term matrix

// int<lower=1,upper=V> w[N]; // word n, this is the COLUMN index of X, LENGTH = N
// int<lower=1,upper=M> doc[N]; // doc ID for word n, this is the ROW index of X, LENGTH = N

vector<lower=0>[K] alpha; // topic prior
vector<lower=0>[V] beta; // word prior
}

parameters {
simplex[K] theta[M]; // topic dist for doc m, MxK matrix
simplex[V] phi[K]; // word dist for topic k, KxV matrix
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
        gamma[k] = log(theta[m, k]) + log(phi[k, v]); // prob(word | theta, phi) = theta * phi 
      }
    target += log_sum_exp(gamma); // likelihood;
  }
}
}
