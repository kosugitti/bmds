data{
  int<lower=2> K; #number of clusters
  int<lower=1> D; #number of dimensions(variables)
  int<lower=1> N; #number of observations
  vector[D] X[N]; #observed data
}

parameters{
  vector[D] mu[K];
  vector<lower=0>[D] sigma[K];
  cov_matrix[D] C[K];
  simplex[K] theta;
}

transformed parameters{
  matrix[D,D] Sigma[K];
  real ps[N,K];

  for(k in 1:K){
    Sigma[k] <- diag_pre_multiply(sigma[k],C[k]);
    Sigma[k] <- diag_post_multiply(Sigma[k],sigma[k]);
  }
  
  for(n in 1:N){
    for(k in 1:K){
      ps[n,k] <- log(theta[k])+multi_normal_log(X[n],mu[k],Sigma[k]);
    }
  }
  
}

model{
  for(n in 1:N){
    increment_log_prob(log_sum_exp(ps[n]));
  }
  for(k in 1:K){
    sigma[k] ~ inv_gamma(0.001,0.001);
  }
}

generated quantities{
  simplex[K] u[N]; #class membership probability
  real log_lik[N];
  for(n in 1:N){
    u[n] <- softmax(to_vector(ps[n]));
    log_lik[n] <- log_sum_exp(ps[n]);
  }
}
