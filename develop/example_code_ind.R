rm(list=ls())
# B-MDS
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanmodel <- stan_model("develop/bmds_ind.stan",model_name="bmds_ind")


library(smacof)
data(helm)
res.helm <- smacofIndDiff(helm, constraint = "indscal")
res.helm$gspace
res.helm$cweights

H <- 16
N <- 10
Npairs <- 45
P <- 2
D <- matrix(unlist(helm),nrow=H,ncol=(N*(N-1)/2),byrow=T)

C <- 4
max.iter <- 2000

### for BMDS
standata <- list(H=H,N=N,Npairs=Npairs,P=P,D=D)
  
fit_vb <- vb(stanmodel,data=standata)
fit_sp <- sampling(stanmodel,data=standata,chain=C,iter=max.iter) #label.switching?

print(fit_sp,digit=2,pars="x")
print(fit_vb,digit=2,pars=c("x"))


#WAIC
library(loo)
log_lik <- rstan::extract(fit_sp)$log_lik
waic(log_lik)
loo(log_lik)

# library(shinystan)
# launch_shinystan(fit_sp)
sp <- rstan::extract(fit_sp,pars="x")$x
### mcmc.list から情報を抜き出す
iter <- length(sp)/(N*P)

#### post MCMC process; For translation-and-reflection issue
#
# case; dim 2
v.flip1 <- matrix(c(1,0,0,1),nrow=2) - (2 * c(1,0) %*% t(c(1,0)))
v.flip2 <- matrix(c(1,0,0,1),nrow=2) - (2 * c(0,1) %*% t(c(0,1)))
V <- list()
idx <- 1
for(v in c(0,pi/2,pi,pi*(3/2))){ #rotation
  Tr <- round(matrix(c(cos(v),sin(v),-sin(v),cos(v)),nrow=2))
  V[[idx]] <- Tr
  idx <- idx + 1
  V[[idx]] <- Tr %*% v.flip1
  idx <- idx + 1
  V[[idx]] <- Tr %*% v.flip2
  idx <- idx + 1
}
V <- unique(V)

Tv <- matrix(nrow=iter,ncol=2)
Tv[,1] <- rep(99,iter)
norm.vec <- c()
maxR <- 50;
for(r in 1:maxR){
  meanx <- apply(sp,c(2,3),mean)
  for(s in 1:iter){
    X <- sp[s,,]
    tmp <- sum(meanx)^2*5 # eps-level init
    for(v in 1:length(V)){
      Xv <- X %*% V[[v]]
      norm.vec[v] <- sum((Xv-meanx)^2)
    }
    sp[s,,] <- sp[s,,] %*% V[[which.min(norm.vec)]]
    Tv[s,2] <- which.min(norm.vec)
  }
  if(sum(Tv[,1]-Tv[,2])==0){
    break;
  }else{
    Tv[,1] <- Tv[,2]
  }
}


################################################# if use vb
# sp <- rstan::extract(fit_vb,pars="x")$x

########### result and plot

configure <- apply(sp,c(2:3),median)
configure.sd <- apply(sp,c(2:3),sd)

plot(configure,type="n")
text(configure,labels=c("Hyogo","Wakayama","Osaka","Nara","Siga","Kyoto"))
x<- seq(-pi,pi,length=100)
for(n in 1:N){
  lines(configure[n,1]+configure.sd[n,1]*cos(x),configure[n,2]+configure.sd[n,2]*sin(x))
}

