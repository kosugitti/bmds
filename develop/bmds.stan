data{
  int<lower=0> N; # data size
  int<lower=0> p; # number of dimensions
  matrix<lower=0>[N,N] D; #disimmirality matrix 
}

parameters{
 vector[N] x[p];
 vector<lower=0>[N*(N-1)/2] phi; #error
 vector<lower=0>[p] lambda ; #Hyper parameter
}

transformed parameters{
  vector<lower=0>[N*(N-1)/2] delta;
  {
    int idx;
    idx <- 1;
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        real di;
        di <- 0.0;
        for(k in 1:p){
          di <- di +(x[i,k]-x[j,k])^2;
        }
        delta[idx] <- sqrt(di);
        idx <- idx + 1;
      }
    }
  }
}

model{
  int idx;
  idx <- 1;
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      D[i,j] ~ normal(delta[idx],phi[idx]);
      phi[idx] ~ gamma(0.001,0.001);
      idx <- idx + 1;
    }
  }
  for(i in 1:p){
    x[i] ~ normal(0,lambda[i]);
    lambda[i] ~ gamma(0.001,0.001);
  }
}
