data {
  int<lower=0> N;  // number of data points
  int<lower=1> D;  // number of dimensions
  int<lower=1> K;  // number of clusters
  vector[D] X[N];  // observations
}

parameters {
  vector[D] mu[K]; // cluster means
  simplex[K] theta;

}
transformed parameters {
  real<upper=0> soft_z[N,K]; // log unnormalized cluster assigns
  for (n in 1:N)
    for (k in 1:K)
      soft_z[n,k] <- log(theta[k]) - 0.5 * dot_self(mu[k] - X[n]);
}
model {
  for (k in 1:K)
    mu[k] ~ normal(0,1);  // prior
  for (n in 1:N)
    increment_log_prob(log_sum_exp(soft_z[n])); // likelihood
}
generated quantities{
  simplex[K] u[N]; #class membership probability
  for(n in 1:N){
    u[n] <- softmax(to_vector(soft_z[n]));
  }
}

