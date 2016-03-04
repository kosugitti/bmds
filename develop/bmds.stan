data{
  int<lower=0> N; # data size
  int<lower=1> Npairs; # number of pairs of data
  int<lower=0> P; # number of dimensions
  row_vector<lower=0>[Npairs] D; #observed dissimirality
}

parameters{
 vector[P] x[N];
 vector<lower=0>[Npairs] phi; #error
 vector<lower=0>[P] invtau2 ; #Hyper parameter
}

transformed parameters{
  vector<lower=0>[Npairs] delta;
  vector<lower=0>[P] tau;  
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
  
  for(p in 1:P)
    tau[p] <- inv_sqrt(invtau2[p]);
}

model{
  for(idx in 1:Npairs){
    D[idx] ~ normal(delta[idx],phi[idx]);
    phi[idx] ~ cauchy(0,2.5);
  }
  for(p in 1:P){
    x[p] ~ normal(0,tau[p]);
  }
  invtau2 ~ gamma(0.001,0.001);
}

generated quantities{
  real log_lik[Npairs];
  for(idx in 1:Npairs)
    log_lik[idx] <- normal_log(D[idx],delta[idx],phi[idx]);
}

