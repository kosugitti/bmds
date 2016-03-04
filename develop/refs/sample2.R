library(bayesmix)
data("fish",package="bayesmix")
x <- fish[,1]
n <- length(x)
K <- 5
m <- 11000
burn <- 1000
model <- BMMmodel(fish,k=K,initialValues=list(S0=2),
                 priors=list(kind='independence',parameter='priorsFish',
                             hierarchical='tau'))
control <- JAGScontrol(variables=c("mu","tau","eta","S"),
                       burn.in=burn,n.iter=m,seed=10)
mcmc <- JAGSrun(fish,model=model,control=control)

J <- 3
mcmc.pars <- array(data=NA,dim=c(m-burn,K,J))
mcmc.pars[,,1]<- mcmc$results[-(1:burn),(n+K+1):(n+2*K)]
mcmc.pars[,,2]<- mcmc$results[-(1:burn),(n+2*K+1):(n+3*K)]
mcmc.pars[,,3]<- mcmc$results[-(1:burn),(n+1):(n+K)]
z <- mcmc$results[-(1:burn),1:n]

library(label.switching)
ls <- label.switching(method="ECR",
                      zpivot=z[mapindex,],
                      z=z,
                      K=K,
                      mcmc=mcmc.pars,
                      constraint=1,
                      data=x)
