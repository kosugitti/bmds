rm(list=ls())
# cmd
iris.dist <-dist(iris[,-5]) #ラベルの部分を除く！
iris.cmd<-cmdscale(iris.dist)
plot(iris.cmd,type='n')
iris.lab <-factor(c(rep("S",50),rep("C",50),rep("V",50)))
text(iris.cmd,labels=iris.lab,col=unclass(iris.lab))

# B-MDS
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanmodel <- stan_model("develop/bmds.stan",model_name="bmds")

citydist <- matrix(c(  0,134, 85,116,118, 60,
                     134,  0, 68, 66,145, 141,
                      85, 68,  0, 32, 83, 75,
                     116, 66, 32,  0, 79, 95,
                     118,145, 83, 79,  0, 63,
                      60,141, 75, 95, 63, 0),nrow=6)
N <- 6
P <- 2
sc <- cmdscale(citydist,P)
# plot(sc)
# text(sc,labels=c("Hyogo","Wakayama","Osaka","Nara","Siga","Kyoto"))

standata <- list(N=N,P=P,D=citydist)
initdata <- list(x=sc)

fit_vb <- vb(stanmodel,data=standata,init=initdata)
fit_sp <- sampling(stanmodel,data=standata,init=list(initdata,initdata,initdata,initdata),chain=4,iter=5000) #label.switching?

fit_sp <- sampling(stanmodel,data=standata,chain=4,iter=5000) #label.switching?

print(fit_sp,digit=2,pars="x")
print(fit_vb,digit=2,pars="x")
# library(shinystan)
# launch_shinystan(fit_sp)
sp <- rstan::extract(fit_sp,pars="x")$x
### mcmc.list から情報を抜き出す
iter <- length(sp)/(N*P)
a <- apply(sp,c(2,3),mean)

sp2 <- as.data.frame(matrix(unlist(sp),nrow=2500,ncol=45))
ret <- matrix(apply(sp2,2,mean)[1:12],ncol=2)
plot(ret,type="n")
text(ret,labels=c("Hyogo","Wakayama","Osaka","Nara","Siga","Kyoto"))

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


