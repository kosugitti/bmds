data{
 int<lower=1> H; #subjects
 int<lower=1> I; #stimuli
 int<lower=1> Ipairs; #pairs of stimuli
 int<lower=1> J; # dimensions
 int<lower=0> K; # latent classes
 row_vector<lower=0>[Ipairs] Y[H]; #observed data
}

parameters{
  simplex[K] theta; #mixing proportions
  matrix[J,K] X[I];
  vector<lower=0>[J] invtau2;
  real<lower=0> invsigma2;
  simplex[J] W0[K,H];
}

transformed parameters{
  matrix<lower=0>[H,Ipairs] d[K];
  vector<lower=0>[J] tau;
  real<lower=0> sigma;
  real ps[H,K];
  row_vector<lower=0>[J] W[K,H];
  for(j in 1:J){
    for(h in 1:H){
      for(k in 1:K){
        W[k,h,j] <- W0[k,h,j];
      }
    }
  }
  for(h in 1:H){
    for(i in 2:I){
      for(ii in 1:(i-1)){
        int t;
        t <- ii + (i-1) * (i-2) /2;
        for(k in 1:K){
          d[k,h,t] <- sqrt(W[k,h,1]*(X[i,1,k]-X[ii,1,k])^2
           + W[k,h,2]*(X[i,2,k]-X[ii,2,k])^2);
          }
      }
    }
  }
  sigma <- inv_sqrt(invsigma2);
  for(j in 1:J){
    tau[j] <- inv_sqrt(invtau2[j]);
  }
  for(h in 1:H){
    for(k in 1:K){
      ps[h,k] <- log(theta[k]+normal_log(Y[h],d[k,h],sigma));
    }
  }
}

model{
  for(i in i:I){
    for(j in 1:J){
      for(k in 1:K){
        X[i,j,k] ~ normal(0,tau[j]);
      }
    }
  }
  invsigma2 ~ gamma(0.001,0.001);
  invtau2 ~ gamma(0.001,0.001);
  for(h in 1:H){
    increment_log_prob(log_sum_exp(ps[h]));
  }
}

generated quantities{
  simplex[K] u[H]; #class membership probability
  matrix[H,Ipairs] Y_pred[K]; #predicted distribution
  for(h in 1:H){
    u[h] <- softmax(to_vector(ps[h]));
  }
  for(k in 1:K){
    for(h in 1:H){
      for(t in 1:Ipairs){
        Y_pred[k,h,t] <- normal_rng(d[k,h,t],sigma);
      }
    }
  }
}
