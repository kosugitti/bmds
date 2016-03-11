data{
  int<lower=1> N; // number of obs
  vector[N] X[N]; //
}

parameters{
  vector[2] xi[N]; // configure
  vector<lower=0>[N] s2[N];
  real<lower=0> mu[N] ; // vM param 1
  real<lower=0> kappa[N]; // vM param 2
}

transformed parameters{
  vector<lower=-3.15,upper=3.15>[N] theta[N];
  vector[N] delta[N];
  vector[N] pi[N];
  vector[N] d[N];

  for(i in 1:N){
    for(j in 1:N){
      theta[i,j] <- atan2(xi[i,2]-xi[j,2],xi[i,1]-xi[j,1]);
      delta[i,j] <- sqrt((xi[i,1]-xi[j,1])^2+(xi[i,2]-xi[j,2])^2);
      pi[i,j] <- exp(von_mises_log(theta[i,j],mu[i],kappa[i]));
      d[i,j] <- (1-pi[i,j])*delta[i,j];
    }
  }
  
}

model{
  for(i in 1:N){
    for(j in 1:N){
      X[i,j] ~ normal(d[i,j],s2[i,j]);
      s2[i,j] ~ cauchy(0,5);
    }
  }
}
