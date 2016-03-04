data{
  int<lower=2> K; #number of clusters
  int<lower=1> N; #number of observations
  real X[N]; #observed data
}

parameters{
  vector[K] mu;
  vector<lower=0>[K] sig2;
  simplex[K] theta;
}

transformed parameters{
  real ps[N,K];
  for(n in 1:N){
    for(k in 1:K){
      ps[n,k] <- log(theta[k])+normal_log(X[n],mu[k],sig2[k]);
    }
  }

}

model{
  for(n in 1:N){
    increment_log_prob(log_sum_exp(ps[n]));
  }
  sig2 ~cauchy(0,2.5);
}
generated quantities{
  simplex[K] u[N]; #class membership probability
  real log_lik[N];
  for(n in 1:N){
    u[n] <- softmax(to_vector(ps[n]));
    log_lik[n] <- log_sum_exp(ps[n]);
  }
  
  
}
