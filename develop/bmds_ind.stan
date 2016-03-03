data{
  int<lower=1> H; # number of subjects
  int<lower=1> Npairs; # number of pairs of data
  int<lower=1> N; # number of stimuli
  int<lower=1> P; # number of dimensions
  row_vector<lower=0>[Npairs] D[H];
}

parameters{
  vector[P] x[N]; # configure
  vector<lower=0>[Npairs] phi; #error
  vector<lower=0>[P] w[H]; # weights
  vector<lower=0>[H] tau; # sd of weights
  vector<lower=0>[P] lambda ; #Hyper parameter
}

transformed parameters{
  vector<lower=0>[Npairs] delta[H];
  for(h in 1:H){
    int idx;
    idx <- 1;
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        real di;
        di <- 0.0;
        for(p in 1:P){
          di <- di + (w[h,p] * (x[i,p]-x[j,p])^2);
        }
      delta[h,idx] <- sqrt(di);
      idx <- idx +1;
      }
    }
  }
}

model{
  for(h in 1:H){
    int idx;
    idx <- 1;
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        D[h,idx] ~ normal(delta[h,idx],phi[idx]);
        idx <- idx + 1;
      }
    }
  }

  for(p in 1:P){
    x[p] ~ normal(0,lambda[p]);
    lambda[p] ~ gamma(0.001,0.001);
  }
  
  for(h in 1:H){
    for(p in 1:P){
      w[h,p] ~ normal(1,tau[h]);
    }
    tau[h] ~ cauchy(0,2.5);
  }
}

generated quantities{
  real log_lik[Npairs];
  for(h in 1:H){
    int idx;
    idx <- 1;
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        log_lik[idx] <- normal_log(D[h,idx],delta[h,idx],phi[idx]);
        idx <- idx +1;
      }
    }
  }
}

