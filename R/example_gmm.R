rm(list=ls())
# B-MDS
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 1 dimensional Mixed Gauss
#  K個の正規分布を混ぜる

K <- 3
mean1 <- 10
mean2 <- -10
mean3 <- 0
sd1 <- 1
sd2 <- 2
sd3 <- 3

N1 <- 500
N2 <- 300
N3 <- 200
N <- N1+N2+N3
X <-c(rnorm(N1,mean1,sd1),rnorm(N2,mean2,sd2),rnorm(N3,mean3,sd3))

### MCMC!
stanmodel <- stan_model("develop/gmm_single.stan",model_name="GMM_Single")
standata <- list(K=K,N=N,X=X)
max.iter <- 3000
burn.in <- 1000
itr <- max.iter - burn.in
C <- 4


# vbは速いしラベルスイッチングとかないから嬉しい
fit_vb <- vb(stanmodel,data=standata)
print(fit_vb,pars=c("mu","sig2","theta"))



# 一つならラベルスイッチングは起きない
fit_sp <- sampling(stanmodel,data=standata,chain=1,iter=max.iter) 
# Rhatも優秀
print(fit_sp,pars=c("mu","sig2","theta"))

# チェインを複数にするとラベルスイッチングが起きる　
fit_sp <- sampling(stanmodel,data=standata,chain=C,iter=max.iter,warmup=burn.in)
# Rhatが悪くなってもくじけちゃいけない
print(fit_sp,pars=c("mu","sig2","theta"))
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
