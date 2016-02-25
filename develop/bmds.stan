data{
  int<lower=0> N; # data size
  int<lower=0> p; # number of dimensions
  matrix<lower=0>[N,N] D; #disimmirality matrix 
}

parameters{
 vector[N] x[p];
 matrix<lower=0>[N,N] phi; #error
 vector<lower=0>[p] lambda ; #Hyper parameter
}

transformed parameters{
  matrix<lower=0>[N,N] delta;
  real di;
  for(i in 1:N){
    for(j in (i+1):N){
      di <- 0;
      for(k in 1:p){
        di <- di +(x[k,i]-x[k,j])^2;
      }
      delta[i,j] <- sqrt(di);
      delta[j,i] <- delta[i,j];
    }
    delta[i,i] <- 0;
  }
}

model{
  for(i in 1:N){
    for(j in (i+1):N){
      D[i,j] ~ normal(delta[i,j],phi[i,j]);
      D[j,i] ~ normal(delta[j,i],phi[j,i]);
      phi[i,j] ~ cauchy(0,5);
      phi[j,i] ~ cauchy(0,5);
    }
  }
  for(i in 1:p){
    x[i] ~ normal(0,lambda[i]);
    lambda[i] ~ cauchy(0,5);
  }
}
