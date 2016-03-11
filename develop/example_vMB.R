rm(list=ls())
# B-MDS
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


example1 <- matrix(c(0,1,5,5,0,1,1,5,0),nrow=3,byrow=T)
example2 <- matrix(c(0,1,1,7,1,0,1,7,7,7,0,1,1,1,7,0),nrow=4)
example3 <- matrix(c(0,7,7,7,1,0,1,7,1,7,0,1,1,1,7,0),nrow=4)
example4 <- matrix(c(0,5,4,5,3,3,2,4,1,1,
                     6,0,4,2,1,2,3,3,4,3,
                     4,4,0,3,3,4,4,5,4,3,
                     4,1,2,0,1,1,4,2,4,3,
                     7,1,2,1,0,1,2,2,2,3,
                     4,3,4,2,3,0,4,4,4,4,
                     4,3,4,4,5,5,0,2,4,2,
                     6,4,4,4,3,4,3,0,4,4,
                     2,3,3,3,3,2,3,3,0,2,
                     4,4,4,5,4,4,4,4,4,0),nrow=10)

dat <- example1/apply(example1,2,max)

ini <- list(xi=cmdscale(dat))

standata <- list(N=ncol(dat),X=dat)
stanmodel <- stan_model("develop/vonMisesBayes.stan",model_name="vonMises")


fit_vb <- vb(stanmodel,data=standata)
fit_sp <- sampling(stanmodel,data=standata)


