data{
  int<lower=0> N; # data size
  int<lower=0> P; # number of dimensions
  matrix<lower=0>[N,N] D; #disimmirality matrix 
}

parameters{
 vector[P] x[N];
 vector<lower=0>[N*(N-1)/2] phi; #error
 vector<lower=0>[P] lambda ; #Hyper parameter
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
        for(p in 1:P){
          di <- di +(x[i,p]-x[j,p])^2;
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
  for(p in 1:P){
    x[p] ~ normal(0,lambda[p]);
    lambda[p] ~ gamma(0.001,0.001);
  }
}

