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
                      60, 14, 75, 95, 63, 0),nrow=6)
N <- 6
sc <- cmdscale(citydist,2)
plot(sc)
text(sc,labels=c("Hyogo","Wakayama","Osaka","Nara","Siga","Kyoto"))

standata <- list(N=N,p=2,D=citydist)

fit_vb <- vb(stanmodel,data=standata)
fit_sp <- sampling(stanmodel,data=standata)
print(fit_vb,digit=2,pars="x")
#print(fit_sp,digit=2)
#library(shinystan)
#launch_shinystan(fit_vb)
sp <- rstan::As.mcmc.list(fit_vb)
sp2 <- as.data.frame(matrix(unlist(sp),nrow=1000,ncol=45))
ret <- matrix(apply(sp2,2,mean)[1:12],ncol=2)
plot(ret)
text(ret,labels=c("Hyogo","Wakayama","Osaka","Nara","Siga","Kyoto"))

