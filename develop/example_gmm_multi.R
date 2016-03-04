rm(list=ls())
# B-MDS
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 3 dimensional Mixed Gauss
data(iris)
dat <- iris[,-5]
K <- 3


### MCMC!
V <- ncol(dat)
N <- nrow(dat)

# stanmodel <- stan_model("develop/refs/soft-k-means.stan",model_name="soft-K-means")
stanmodel <- stan_model("develop/gmm_softKmeans_simplex.stan",model_name="softKmeans_simplex")
# stanmodel <- stan_model("develop/gmm_softKmeans_multi.stan",model_name="softKmeans_multi") # divergent?!
# stanmodel <- stan_model("develop/gmm_simplex_multi.stan",model_name="simplex_multi") # Too Slow!
max.iter <- 3000
burn.in <- 1000
itr <- max.iter - burn.in
C <- 4
standata <- list(K=K,N=N,D=V,X=dat)

set.seed(42)
init.kmeans <- kmeans(dat,K,algorithm = "MacQueen")
kmeans.center <- init.kmeans$centers
covs <- list()
for(i in 1:K){
  sub <- subset(dat,init.kmeans$cluster==i)
  covs[[i]] <- cov(sub)
}
init <- list(mu=kmeans.center,Sigma=covs)

# vbは速いしラベルスイッチングとかないから嬉しい
fit_vb <- vb(stanmodel,data=standata)
print(fit_vb,pars="mu")



# 一つならラベルスイッチングは起きない
fit_sp <- sampling(stanmodel,data=standata,init=list(init),chain=1,iter=max.iter)
# Rhatも優秀
print(fit_sp,pars=c("mu"))

# チェインを複数にするとラベルスイッチングが起きる　
init.list <- rep(init,C)
fit_sp <- sampling(stanmodel,data=standata,chain=C,iter=max.iter,warmup=burn.in)
# Rhatが悪くなってもくじけちゃいけない
print(fit_sp,pars=c("mu"))
# 目で見ると綺麗なラベルスイッチングが確認できるよ。
traceplot(fit_sp,pars="mu")


################################################################### label.switching関数の準備です
allocK <- rstan::extract(fit_sp,pars="u")$u
zChain <- matrix(ncol=N,nrow=itr*C)
for(i in 1:(itr*C)){
  zChain[i,] <- apply(allocK[i,,],1,which.max)
}


mcmc.params <- rstan::extract(fit_sp,pars=c("mu","sig2","theta"))

# mcmc * num.of.clusters * num.of.params
mcmc.pars <- array(data=NA,dim=c(itr*C,K,3))
mcmc.pars[,,1]<- mcmc.params$mu
mcmc.pars[,,2]<- mcmc.params$sig2
mcmc.pars[,,3]<- mcmc.params$theta

# この関数はSJW法を使う時に必要になってくる。completeオプションに関数を渡さないといけないので！
complete.normal.loglikelihood<-function(x,z,pars){
  #	x: data (size = n)
  #	z: allocation vector (size = n)
  #	pars: K\times J vector of normal mixture parameters:
  #		pars[k,1] = mean of the k-normal component
  #		pars[k,2] = variance of the k-normal component
  #		pars[k,3] = weight of the k-normal component
  #			k = 1,...,K
  g <- dim(pars)[1] #K (number of mixture components)
  n <- length(x)	#this denotes the sample size
  logl<- rep(0, n)	
  logpi <- log(pars[,3])
  mean <- pars[,1]
  sigma <- sqrt(pars[,2])
  logl<-logpi[z] + dnorm(x,mean = mean[z],sd = sigma[z],log = TRUE)
  return(sum(logl))
}


# パッケージのお出まし
library(label.switching)
# 全方法試す
set <- c("STEPHENS", "PRA", "ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2", "AIC", "SJW","DATA-BASED")

#ピボットは対数尤度が一番小さいやつにしてみた（なんでもいいと思う）
log_liks <- rstan::extract(fit_sp,pars="log_lik")$log_lik
sort(table(apply(log_liks,2,which.min)),decreasing = T)
pivot <- 158

# さあ実行
ls <- label.switching(method = set,
                      zpivot = zChain[pivot,],
                      prapivot = mcmc.pars[pivot, , ],                      
                      z = zChain,
                      K = 3, 
                      p = allocK,
                      mcmc=mcmc.pars,
                      complete=complete.normal.loglikelihood,
                      data = X)
# 各推定法で結果が一致したかどうか
ls$similarity
# 推定された所属クラスタ番号
ls$clusters
