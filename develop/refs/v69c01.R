# The following packages are required
library("label.switching")  #(version 1.3)       
library("bayesmix")  # rjags is required here, see: http://ifas.jku.at/gruen/BayesMix    
library("RColorBrewer")
library("mvtnorm")
dir.create("labelSwitchingFigures")  #  all output (pdf-plots) will be saved here  
setwd("labelSwitchingFigures/")

# *** Note that using source('v69c03.R') takes ~ 2.5-3 hours (Linux machine)

# ===============================================================================================
# (1) Fishery data example
# ===============================================================================================

data("fish", package = "bayesmix")
x <- fish[, 1]
n <- length(x)  # # sample size: n = 256
K <- 5  #  # number of components
# produce an MCMC sample using bayesmix package
m <- 11000  #  # define number of iterations
burn <- 1000  # # define burn-in period
model <- BMMmodel(fish, k = K, initialValues = list(S0 = 2), priors = list(kind = "independence", 
  parameter = "priorsFish", hierarchical = "tau"))
control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"), burn.in = burn, n.iter = m, 
  seed = 10)
mcmc <- JAGSrun(fish, model = model, control = control)
# write output to a convenient format for label.switching package: MCMC
# parameters should be saved as m\times K\times J array
J <- 3  #   # three different types of parameters:
# j=1 corresponds to simulated means
# j=2 corresponds to simulated variances
# j=3 corresponds to simulated weights

mcmc.pars <- array(data = NA, dim = c(m - burn, K, J))
mcmc.pars[, , 1] <- mcmc$results[-(1:burn), (n + K + 1):(n + 2 * K)]
mcmc.pars[, , 2] <- mcmc$results[-(1:burn), (n + 2 * K + 1):(n + 3 * K)]
mcmc.pars[, , 3] <- mcmc$results[-(1:burn), (n + 1):(n + K)]
# save the generated allocation variables to m\times n array z:
z <- mcmc$results[-(1:burn), 1:n]
colnames(z) <- NULL
m <- m - burn
# Compute necessary information that will be additional input to the
# label.switching package.  Define the complete log-likelihood function of the
# normal mixture
complete.normal.loglikelihood <- function(x, z, pars) {
  g <- dim(pars)[1]
  n <- length(x)
  logl <- rep(0, n)
  logpi <- log(pars[, 3])
  mean <- pars[, 1]
  sigma <- sqrt(pars[, 2])
  logl <- logpi[z] + dnorm(x, mean = mean[z], sd = sigma[z], log = T)
  return(sum(logl))
}
# computing complete loglikelihoods in order to find pivot
iter <- 1
mapindex <- 1
x <- fish[, 1]
pars <- mcmc.pars[iter, , ]
maxL <- complete.normal.loglikelihood(x, z[iter, ], pars)
for (iter in 2:m) {
  pars <- mcmc.pars[iter, , ]
  logL <- complete.normal.loglikelihood(x, z[iter, ], pars)
  if (logL > maxL) {
    maxL <- logL
    mapindex <- iter
  }
}
print(paste("complete likelihood pivot = ", mapindex))
mapindex <- mapindex
zmap <- z[mapindex, ]
# computing allocation probabilities for stephens method and ECR-ITERATIVE-2
p <- array(data = NA, dim = c(m, n, K))
for (iter in 1:m) {
  for (i in 1:n) {
    kdist <- mcmc.pars[iter, , 3] * dnorm(x[i], mcmc.pars[iter, , 1], sqrt(mcmc.pars[iter, 
      , 2]))
    skdist <- sum(kdist)
    for (j in 1:K) {
      p[iter, i, j] <- kdist[j]/skdist
    }
  }
}

# Run label.switching command using all methods.

# The pivot for default ECR algorithm will be the allocation `mapindex` that
# corresponds to the maximum of complete logL estimate: zpivot=z[mapindex,].

# The pivot for PRA algorithm will be the parameters that correspond to the same
# iteration: prapivot = mcmc.pars[mapindex,,].

# The SJW method will be initialized using this iteration as well: sjwinit =
# mapindex.  The complete log-likelihood is defined as: complete =
# complete.normal.loglikelihood

set <- c("STEPHENS", "PRA", "ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2", "SJW", "AIC", 
  "DATA-BASED")
ls <- label.switching(method = set, zpivot = zmap, z = z, K = K, prapivot = mcmc.pars[mapindex, 
  , ], p = p, complete = complete.normal.loglikelihood, mcmc.pars, data = x, sjwinit = mapindex)

ls$similarity  #similarity of single best clusterings:      
ls$timings  #timings for each method

t <- 5
thinning <- t + seq(1, m, by = t) - 1  # #select every t-th iteration for plotting the mcmc 
colors <- brewer.pal(K, name = "Set1")  # #define color pallete     
pal <- colorRampPalette(colors)

# plot the simulated means (label switching phenomenon)
image.width <- 6
image.height <- 6
pointsize <- 20
mymargin <- c(4, 4, 2, 0.5)
pdf(file = "fish-raw-means.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(mcmc.pars[thinning, , 1], type = "l", lty = 1, col = pal(K), ylab = "means", 
  xlab = "iteration", main = "Raw MCMC output")
dev.off()

# Reordered outputs
pdf(file = "fish-ecr.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$ECR)$output[thinning, , 1], type = "l", 
  xlab = "iteration", main = "ECR", ylab = "means", lty = 1, col = pal(K))
dev.off()
pdf(file = "fish-ecr-iterative-2.pdf", width = image.width, height = image.height, 
  pointsize = pointsize)
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$"ECR-ITERATIVE-2")$output[thinning, 
  , 1], type = "l", xlab = "iteration", main = "ECR-iterative-2", ylab = "means", 
  lty = 1, col = pal(K))
dev.off()
pdf(file = "fish-ecr-iterative-1.pdf", width = image.width, height = image.height, 
  pointsize = pointsize)
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$"ECR-ITERATIVE-1")$output[thinning, 
  , 1], type = "l", xlab = "iteration", main = "ECR-iterative-1", ylab = "means", 
  lty = 1, col = pal(K))
dev.off()
pdf(file = "fish-pra.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$PRA)$output[thinning, , 1], type = "l", 
  xlab = "iteration", main = "PRA", ylab = "means", lty = 1, col = pal(K))
dev.off()
pdf(file = "fish-stephens.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$STEPHENS)$output[thinning, , 1], 
  type = "l", xlab = "iteration", main = "STEPHENS", ylab = "means", lty = 1, col = pal(K))
dev.off()
pdf(file = "fish-sjw.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$SJW)$output[thinning, , 1], type = "l", 
  xlab = "iteration", main = "SJW", ylab = "means", lty = 1, col = pal(K))
dev.off()
pdf(file = "fish-aic.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$AIC)$output[thinning, , 1], type = "l", 
  xlab = "iteration", main = "AIC", ylab = "means", lty = 1, col = pal(K))
dev.off()
pdf(file = "fish-dataBased.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$"DATA-BASED")$output[thinning, , 
  1], type = "l", xlab = "iteration", main = "DATA-BASED", ylab = "means", lty = 1, 
  col = pal(K))
dev.off()

# histogram of the observed data
pdf(file = "fish-hist.pdf", width = image.height, height = image.height)
par(mar = c(4, 4, 2, 0))
hist(x, 50, xlab = "length", main = "Fishery data")
dev.off()

# ===============================================================================================
# (2) Fetal lamb data example
# ===============================================================================================

data("fetal-lamb", package = "label.switching")
x <- lamb
K <- 4  #  #set number of components
n <- length(x)  # #sample size

# GIBBS SAMPLER
set.seed(9583803)
# set number of iterations and burn in
m <- 20000
burn <- 10000
# fruhwirth--schnatter prior on the state transition matrix
dir.prior <- array(data = 1/(K - 1), dim = c(K, K))
diag(dir.prior) <- rep(2, K)
log.dir.prior <- function(prob) {
  logprior <- 0
  for (k in 1:K) {
    logprior <- logprior + sum((dir.prior[k, ] - 1) * log(prob[k, ]))
  }
  return(logprior)
}
# fruhwirth-schnatter gamma prior parameters for lambda_k, k=1,...,K
alpha <- 1
beta <- 0.5
psim <- array(data = NA, dim = c(m, K, K))
theta <- array(data = NA, dim = c(m, K))
z <- array(data = NA, dim = c(m, n))
eigenvectors <- array(data = NA, dim = c(m, K))
# initialization
iter <- 1
theta[iter, ] <- rgamma(K, shape = alpha, rate = beta)
for (k in 1:K) {
  psim[iter, k, ] <- rgamma(K, shape = dir.prior[k, ], rate = 1)
  psim[iter, k, ] <- psim[iter, k, ]/sum(psim[iter, k, ])
}
eigenvectors[iter, ] <- Re(eigen(t(psim[iter, , ]))$vectors[, 1]/sum(eigen(t(psim[iter, 
  , ]))$vectors[, 1]))
z[iter, 1] <- sample(1:K, 1, prob = eigenvectors[iter, ])
for (i in 2:n) {
  z[iter, i] <- sample(1:K, 1, prob = psim[iter, z[iter, i - 1], ])
}
post <- sum(dpois(x, theta[iter, z[iter, ]], log = T))
for (i in 2:n) {
  post <- post + log(psim[iter, z[iter, i - 1], z[iter, i]])
}
post <- post + log(eigenvectors[iter, z[iter, 1]]) + sum(dgamma(theta[iter, ], shape = alpha, 
  rate = beta, log = T)) + log.dir.prior(psim[iter, , ])
maxval <- post
mapindex <- iter
post.non <- 0
for (i in 1:n) {
  post.non <- post.non + log(sum(eigenvectors[iter, ] * dpois(x[i], theta[iter, 
    ])))
}
maxval.non <- post.non
mapindex.non <- iter
for (iter in 2:m) {
  nij <- array(data = 0, dim = c(K, K))
  nj <- sj <- numeric(K)
  # update z
  i <- 1
  # compute probs
  prob <- eigenvectors[iter - 1, ] * dpois(x[i], theta[iter - 1, ]) * psim[iter - 
    1, , z[iter - 1, i + 1]]
  prob <- prob/sum(prob)
  z[iter, i] <- sample(1:K, 1, prob = prob)
  nj[z[iter, i]] <- nj[z[iter, i]] + 1
  sj[z[iter, i]] <- sj[z[iter, i]] + x[i]
  for (i in 2:(n - 1)) {
    prob <- psim[iter - 1, z[iter, i - 1], ] * dpois(x[i], theta[iter - 1, ]) * 
      psim[iter - 1, , z[iter - 1, i + 1]]
    prob <- prob/sum(prob)
    z[iter, i] <- sample(1:K, 1, prob = prob)
    nj[z[iter, i]] <- nj[z[iter, i]] + 1
    sj[z[iter, i]] <- sj[z[iter, i]] + x[i]
    nij[z[iter, i - 1], z[iter, i]] <- nij[z[iter, i - 1], z[iter, i]] + 1
  }
  i <- n
  prob <- psim[iter - 1, z[iter, i - 1], ] * dpois(x[i], theta[iter - 1, ])
  prob <- prob/sum(prob)
  z[iter, i] <- sample(1:K, 1, prob = prob)
  nj[z[iter, i]] <- nj[z[iter, i]] + 1
  sj[z[iter, i]] <- sj[z[iter, i]] + x[i]
  nij[z[iter, i - 1], z[iter, i]] <- nij[z[iter, i - 1], z[iter, i]] + 1
  # update matrix of probabilities
  for (k in 1:K) {
    psim[iter, k, ] <- rgamma(K, shape = dir.prior[k, ] + nij[k, ], rate = 1)
    psim[iter, k, ] <- psim[iter, k, ]/sum(psim[iter, k, ])
  }
  # compute eigenvector
  eigenvectors[iter, ] <- Re(eigen(t(psim[iter, , ]))$vectors[, 1]/sum(eigen(t(psim[iter, 
    , ]))$vectors[, 1]))
  # update Poisson parameters
  theta[iter, ] <- rgamma(K, shape = alpha + sj, rate = beta + nj)
  # compute complete posterior density
  post <- sum(dpois(x, theta[iter, z[iter, ]], log = T))
  for (i in 2:n) {
    post <- post + log(psim[iter, z[iter, i - 1], z[iter, i]])
  }
  post <- post + log(eigenvectors[iter, z[iter, 1]]) + sum(dgamma(theta[iter, ], 
    shape = alpha, rate = beta, log = T)) + log.dir.prior(psim[iter, , ])
  if ((post > maxval) & (iter > burn)) {
    maxval <- post
    mapindex <- iter
  }
  
  # compute observed loglikelihood (in case that additional pivots are required)
  # post.non<-0
  # for(i in 1:n){
  # post.non<-post.non+log(sum(eigenvectors[iter,]*dpois(x[i],theta[iter,])))
  # }
  # if((post.non>maxval.non)&(iter>burn)){maxval.non<-post.non;mapindex.non<-iter}
  # if(iter%%500==0){matplot(theta[1:iter,],type='l')}
  if (iter%%500 == 0) {
    matplot(log(theta[1:iter, ]), type = "l", xlab = "log-Mean")
  }
}

# END OF GIBBS SAMPLER
# Rearrange MCMC output to (m,K,J) array
# pars: dim KxJ
# J = 1 + K + 1
# j = 1 corresponds to theta
# j = 2,...,K + 1 corresponds to p[,j]
# p_{ij} = pars[i,j+1]
# j = K + 2 corresponds to eigenvectors
J <- 1 + K + 1
mcmc.pars <- array(data = NA, dim = c(m - burn, K, J))
mcmc.pars[, , 1] <- theta[-(1:burn), ]  # # poisson means      
for (k in 1:K) {
  mcmc.pars[, , 1 + k] <- psim[-(1:burn), , k]
}
mcmc.pars[, , K + 2] <- eigenvectors[-(1:burn), ]
z <- z[-(1:burn), ]
mapindex <- mapindex - burn
# mapindex.non<-mapindex.non -burn

# define complete log-likelihood function for SJW algorithm

complete.hmm.poisson.loglikelihood <- function(x, z, pars) {
  post <- sum(dpois(x, pars[z, 1], log = T))
  logprobs <- log(pars[, 2:(K + 1)])
  ev <- pars[, K + 2]
  for (i in 2:n) {
    post <- post + logprobs[z[i - 1], z[i]]
  }
  post <- post + log(ev[z[1]])
  return(post)
}

# finding the P matrix for Stephens algorithm
p <- array(data = NA, dim = c(m - burn, n, K))
for (iter in (burn + 1):m) {
  for (i in 1:n) {
    kdist <- eigenvectors[iter, ] * dpois(x[i], theta[iter, ])
    skdist <- sum(kdist)
    for (j in 1:K) {
      p[iter - burn, i, j] <- kdist[j]/skdist
    }
  }
}

# Apply relabelling algorithms
set <- c("STEPHENS", "PRA", "ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2", "SJW", "AIC", 
  "DATA-BASED")
ls <- label.switching(method = set, zpivot = z[mapindex, ], z = z, K = K, prapivot = mcmc.pars[mapindex, 
  , ], p = p, complete = complete.hmm.poisson.loglikelihood, mcmc = mcmc.pars, 
  data = x)

# retrieve the number of observations assigned to each cluster:
frequency <- apply(ls$clusters, 1, function(y) {
  freq <- numeric(K)
  for (j in 1:K) {
    freq[j] <- length(which(y == j))
  }
  return(freq)
})
rownames(frequency) <- 1:K
frequency

m <- m - burn
t <- 5
thinning <- t + seq(1, m, by = t) - 1

colors <- brewer.pal(K, name = "Set1")
pal <- colorRampPalette(colors)

ylimit <- c(-6, 2)
image.width <- 6
image.height <- 6
pointsize <- 20
mymargin <- c(4, 4, 2, 0.5)

# plot the raw output of poisson log-means (label switching)
pdf(file = "lamb-raw-means.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(log(mcmc.pars[thinning, , 1]), type = "l", lty = 1, ylab = "log-means", xlab = "iteration", 
  main = "Raw MCMC output", col = pal(K), ylim = ylimit)
dev.off()
# plot the reordered output of poisson means according to each method
pdf(file = "lamb-ecr.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(log(permute.mcmc(mcmc.pars, ls$permutations$ECR)$output[thinning, , 1]), 
  type = "l", xlab = "iteration", main = "ECR", ylab = "log-means", lty = 1, col = pal(K), 
  ylim = ylimit)
dev.off()
pdf(file = "lamb-ecr-iterative-1.pdf", width = image.width, height = image.height, 
  pointsize = pointsize)
par(mar = mymargin)
matplot(log(permute.mcmc(mcmc.pars, ls$permutations$"ECR-ITERATIVE-1")$output[thinning, 
  , 1]), type = "l", xlab = "iteration", main = "ECR-iterative-1", ylab = "log-means", 
  lty = 1, col = pal(K), ylim = ylimit)
dev.off()
pdf(file = "lamb-ecr-iterative-2.pdf", width = image.width, height = image.height, 
  pointsize = pointsize)
par(mar = mymargin)
matplot(log(permute.mcmc(mcmc.pars, ls$permutations$"ECR-ITERATIVE-2")$output[thinning, 
  , 1]), type = "l", xlab = "iteration", main = "ECR-iterative-2", ylab = "log-means", 
  lty = 1, col = pal(K), , ylim = ylimit)
dev.off()
pdf(file = "lamb-pra.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(log(permute.mcmc(mcmc.pars, ls$permutations$PRA)$output[thinning, , 1]), 
  type = "l", xlab = "iteration", main = "PRA", ylab = "log-means", lty = 1, col = pal(K), 
  ylim = ylimit)
dev.off()
pdf(file = "lamb-stephens.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(log(permute.mcmc(mcmc.pars, ls$permutations$STEPHENS)$output[thinning, , 
  1]), type = "l", xlab = "iteration", main = "STEPHENS", ylab = "log-means", lty = 1, 
  col = pal(K), ylim = ylimit)
dev.off()
pdf(file = "lamb-sjw.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(log(permute.mcmc(mcmc.pars, ls$permutations$SJW)$output[thinning, , 1]), 
  type = "l", xlab = "iteration", main = "SJW", ylab = "log-means", lty = 1, col = pal(K), 
  ylim = ylimit)
dev.off()
pdf(file = "lamb-aic.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(log(permute.mcmc(mcmc.pars, ls$permutations$AIC)$output[thinning, , 1]), 
  type = "l", xlab = "iteration", main = "AIC", ylab = "log-means", lty = 1, col = pal(K), 
  ylim = ylimit)
dev.off()
pdf(file = "lamb-data-based.pdf", width = image.width, height = image.height, pointsize = pointsize)
par(mar = mymargin)
matplot(log(permute.mcmc(mcmc.pars, ls$permutations$"DATA-BASED")$output[thinning, 
  , 1]), type = "l", xlab = "iteration", main = "DATA-BASED", ylab = "log-means", 
  lty = 1, col = pal(K), ylim = ylimit)
dev.off()

# time series of the observed data
pdf(file = "lamb-time.pdf", width = image.height, height = image.height)
par(mar = c(4, 4, 2, 0))
plot(x, type = "l", ylab = "number of movements", xlab = "time", main = "Fetal lamb data", 
  bty = "l")
dev.off()

# ===============================================================================================
# (3) Mixture of bivariate normals, simulated data
# ===============================================================================================

# defining necessary functions:

# gibbsSampler: Gibbs Sampling function

# complete.bivariate.normal.loglikelihood: complete log-likelihood for the
# bivariate case

# normDens: log-pdf function for the mixture

# complete log-likelihood for the bivariate mixture
complete.bivariate.normal.loglikelihood <- function(x, z, pars) {  
  # x: denotes the n data points
  # z: denotes an allocation vector (size=n)
  # pars: K\times 3 vector of means,variance, weights
  # pars[k,1]: corresponds to the mean of component k
  # pars[k,2]: corresponds to the variance of component k
  # pars[k,3]: corresponds to the weight of component k
  g <- dim(pars)[1]
  n <- dim(x)[1]
  logl <- rep(0, n)
  logpi <- log(pars[, 6])
  my_mean <- pars[, 1:2]
  d <- 2
  cov.mat <- array(data = NA, dim = c(d, d))
  logl <- logpi[z]
  for (k in 1:g) {
    ind <- which(z == k)
    if (length(ind) > 0) {
      diag(cov.mat) <- pars[k, 3:4]
      cov.mat[1, 2] <- cov.mat[2, 1] <- pars[k, 5]
      logl[ind] <- logl[ind] + dmvnorm(x[ind, ], mean = my_mean[k, ], sigma = cov.mat, 
        log = TRUE)
    }
  }
  return(sum(logl))
}

# log-pdf function for the bivariate mixture
normDens <- function(x, pars) {
  g <- dim(pars)[1]
  pi <- pars[, 6]
  my_mean <- pars[, 1:2]
  d <- 2
  cov.mat <- array(data = NA, dim = c(d, d))
  n <- dim(x)[1]
  if (is.null(n)) {
    n <- 1
  }
  f <- rep(0, n)
  for (k in 1:g) {
    diag(cov.mat) <- pars[k, 3:4]
    cov.mat[1, 2] <- cov.mat[2, 1] <- pars[k, 5]
    f <- f + pi[k] * dmvnorm(x, mean = my_mean[k, ], sigma = cov.mat)
  }
  return(log(f))
}

# Gibbs Sampler
gibbsSampler <- function(iterations, K, x, burn) {
  colors <- brewer.pal(K, name = "Set1")  # #pretty only when K < 10    
  pal <- colorRampPalette(colors)
  d <- dim(x)[2]
  burn <- floor(abs(burn))
  iterations <- floor(abs(iterations))
  if (d != 2) {
    stop(cat(paste("parts of the function are applicable only to bivariate data"), 
      "\n"))
  }
  if (burn > iterations) {
    stop(cat(paste("burn-in period is not valid"), "\n"))
  }
  if (K != floor(K)) {
    stop(cat(paste("K should be integer > 1"), "\n"))
  }
  if (K < 2) {
    stop(cat(paste("K should be larger than 1"), "\n"))
  }
  if (dim(x)[1] < 2) {
    stop(cat(paste("number of observations should be > 2"), "\n"))
  }
  if (iterations < 2) {
    stop(cat(paste("number of iterations should be > 1"), "\n"))
  }
  # Prior distributions, precision matrix for Wishart prior: t_k ~
  # W(priorSigma_k,v)
  v <- rep(d + 1, K)
  invPriorSigma <- cov.mat <- priorSigma <- array(data = 0, dim = c(d, d, K))
  for (k in 1:K) {
    for (j in 1:d) {
      priorSigma[j, j, k] <- 0.5
    }
    invPriorSigma[, , k] <- solve(priorSigma[, , k])
  }
  # Normal distribution hyperparamaters: \mu_k|Sigma ~ N (m_k,b*Sigma)
  b <- rep(0.05, K)
  m <- array(data = NA, dim = c(d, K))
  for (j in 1:d) {
    m[j, ] <- mean(x)
  }
  # Dirichlet prior parameters
  a <- rep(1, K)
  p <- array(data = NA, dim = c(iterations, K))  #  #mixture weights
  mu <- array(data = NA, dim = c(iterations, d, K))  # #mixture means     
  t <- array(data = NA, dim = c(iterations, d, d, K))  # #mixture inverse variance  
  zChain <- array(data = NA, dim = c(iterations, n))  # #mixture allocations   
  allocProbs <- array(data = 0, dim = c(iterations, n, K))  # #classification probabilities  
  # Initialization:
  iter <- 1
  for (k in 1:K) {
    h <- array(priorSigma[, , k], dim = c(d, d))
    t[iter, , , k] <- rWishart(1, df = v[k], h)
    cov.mat[, , k] <- solve(t[iter, , , k])
    h <- array(cov.mat[, , k], dim = c(d, d))
    mu[iter, , k] <- rmvnorm(1, mean = m[, k], sigma = h/b[k])
  }
  p[iter, ] <- rgamma(K, shape = a, rate = 1)
  p[iter, ] <- p[iter, ]/sum(p[iter, ])
  par(mfrow = c(1, d))
  z <- numeric(n)
  probs <- numeric(K)
  zChain[1, ] <- rep(1, n)
  allocProbs[1, , ] <- array(1, dim = c(n, K))
  # Updates
  par(mfrow = c(1, d))
  for (iter in 2:iterations) {
    s1 <- numeric(K)
    s2 <- array(data = 0, dim = c(d, K))
    s3 <- invPriorSigma
    # update z
    for (i in 1:n) {
      for (k in 1:K) {
        h <- array(cov.mat[, , k], dim = c(d, d))
        probs[k] <- p[iter - 1, k] * dmvnorm(x[i, ], mean = mu[iter - 1, 
          , k], sigma = h)
      }
      allocProbs[iter, i, ] <- probs
      z[i] <- sample(1:K, 1, prob = probs)
      s1[z[i]] <- s1[z[i]] + 1
      s2[, z[i]] <- s2[, z[i]] + x[i, ]
      s3[, , z[i]] <- s3[, , z[i]] + x[i, ] %*% t(x[i, ])
    }
    # update mu and t
    for (k in 1:K) {
      if (s1[k] > 0) {
        s3[, , k] <- s3[, , k] - (s2[, k] %*% t(s2[, k]))/s1[k] + (b[k] * 
          s1[k]/(b[k] + s1[k])) * (s2[, k]/s1[k] - m[, k]) %*% t(s2[, k]/s1[k] - 
          m[, k])
      }
      wish.matrix <- solve(s3[, , k])
      wish.matrix <- array(wish.matrix, dim = c(d, d))
      t[iter, , , k] <- rWishart(1, df = v[k] + s1[k], wish.matrix)
      cov.mat[, , k] <- solve(t[iter, , , k])
      h <- array(cov.mat[, , k], dim = c(d, d))
      mm <- (b[k] * m[, k] + s2[, k])/(b[k] + s1[k])
      mu[iter, , k] <- rmvnorm(1, mean = mm, sigma = h/(b[k] + s1[k]))
    }
    # update p
    p[iter, ] <- rgamma(K, shape = a + s1, rate = 1)
    p[iter, ] <- p[iter, ]/sum(p[iter, ])
    zChain[iter, ] <- z
    if (iter%%100 == 0) {
      for (i in 1:d) {
        matplot(mu[1:iter, i, ], type = "l", ylim = range(x), col = pal(K))
        abline(h = mu.real[i, ])
      }
    }
  }
  iterations <- iterations - burn
  p <- p[-(1:burn), ]
  zChain <- zChain[-(1:burn), ]
  t <- t[-(1:burn), , , ]
  mu <- mu[-(1:burn), , ]
  allocProbs <- allocProbs[-(1:burn), , ]
  # invert t to obtain covariance matrices
  inv.t <- t
  for (iter in 1:iterations) {
    for (k in 1:K) {
      inv.t[iter, , , k] <- solve(t[iter, , , k])
    }
  }
  t <- 0
  # compute likelihood values for pivot and normalize classification probs (this
  # could be done inside the mcmc updates as well.)
  lValues <- numeric(iterations)
  mcmc.pars <- array(data = NA, dim = c(iterations, K, 6))
  for (iter in 1:iterations) {
    
    # Note that we rearranged the MCMC output to (iterations,K,J) array
    # mcmc: iterations x K x J
    # J = d + d*(d+1)/2 + 1 is the number of different parameter types
    # j = 1 corresponds to mu[,1] (mean of x[,1])
    # j = 2 corresponds to mu[,2] (mean of x[,2])
    # j = 3 corresponds to Sigma[,1,1] (variance of x[,1])
    # j = 4 corresponds to Sigma[,2,2] (variance of x[,2])
    # j = 5 corresponds to Sigma[,1,2] (covariance of (x[,1],x[,2]))
    # j = 6 corresponds to mixture weights
    mcmc.pars[iter, , ] <- cbind(t(mu[iter, , ]), inv.t[iter, 1, 1, ], inv.t[iter, 
      2, 2, ], inv.t[iter, 1, 2, ], p[iter, ])
    lValues[iter] <- complete.bivariate.normal.loglikelihood(x, zChain[iter, 
      ], mcmc.pars[iter, , ])
    nc <- apply(allocProbs[iter, , ], 1, sum)
    allocProbs[iter, , ] <- allocProbs[iter, , ]/nc
  }
  mapindex <- order(lValues, decreasing = TRUE)[1]
  results <- list(mcmc.pars, mapindex, zChain, allocProbs)
  names(results) <- c("mcmc", "MLindex", "z", "p")
  return(results)
}

# Example with simulated dataset 1 (K = 4)

# Simulation set-up for true parameter values
# K: number of components
# n: number of observations
# d: dimensionality of multivariate normal distribution
# x: data
# mu.real: real values for the means
# t.real: real values for inverse covariance matrices
# z.real: real allocation of observations among the K mixture components
# p.real: real values for the weights of the mixture
set.seed(999)
K <- 4
n <- 100
d <- 2
x <- array(data = NA, dim = c(n, d))
p.real <- rep(1, K)
myDistance <- 2.5
mu.real <- array(data = NA, dim = c(d, K))
tr <- seq(0, 2 * pi, length = K + 1)
tr <- tr[-(K + 1)]
mu.real[1, ] <- cos(tr)
mu.real[2, ] <- sin(tr)
mu.real <- myDistance * mu.real
t.real <- array(data = 0, dim = c(d, d, K))
for (k in 1:K) {
  for (j in 1:d) t.real[j, j, k] <- 1
}
if (d > 1) {
  t.real[1, 2, 1] <- 0
  t.real[2, 1, 1] <- 0
}
z.real <- sample(1:K, n, prob = p.real, replace = TRUE)
x <- array(data = NA, dim = c(n, d))
for (k in 1:K) {
  index <- which(z.real == k)
  diaspora <- solve(t.real[, , k])
  x[index, ] <- rmvnorm(length(index), mean = mu.real[, k], sigma = diaspora, method = c("eigen", 
    "svd", "chol"), pre0.9_9994 = FALSE)
}
plot(x, col = z.real)  #scatterplot of simulated data with true cluster labels
# Run Gibbs Sampler
iterations <- 11000
burn <- 1000
gs <- gibbsSampler(iterations = iterations, K = 4, x = x, burn = burn)  #run-time: ~25minutes
iterations <- iterations - burn
zChain <- gs$z
mcmc.pars <- gs$mcmc
pivot <- gs$MLindex
allocProbs <- gs$p
gs <- 0

# apply user-defined constraint to: mu_1 - 2mu_2
newMCMC <- array(data = NA, dim = c(iterations, K, 7))
newMCMC[, , 1:6] <- mcmc.pars
for (k in 1:K) {
  newMCMC[, k, 7] <- mcmc.pars[, k, 1] - 2 * mcmc.pars[, k, 2]
}
newConstraint <- aic(newMCMC, constraint = 7)
# run label.switching.  In order to compare all methods with true clusters:
# groundTruth = z.real
set <- c("STEPHENS", "PRA", "ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2", "SJW", "AIC", 
  "DATA-BASED", "USER-PERM")
ls <- label.switching(method = set, zpivot = zChain[pivot, ], z = zChain, K = K, 
  p = allocProbs, prapivot = mcmc.pars[pivot, , ], complete = complete.bivariate.normal.loglikelihood, 
  mcmc = mcmc.pars, data = x, sjwinit = pivot, groundTruth = z.real, userPerm = newConstraint$permutations)

# get the similarity of single best clusterings (and check the similarity of
# estimate with the ground truth and among them)
ls$similarity

colors <- brewer.pal(K, name = "Set1")
pal <- colorRampPalette(colors)
# random sample of 500 iterations for clear plotting
thinning <- sample(1:iterations, size = 500, replace = FALSE)
image.width <- 6
image.height <- 6
pointsize <- 20
mymargin <- c(4, 4, 2, 0.5)

# raw MCMC output for the means (with label switching)
pdf(file = "bivariate-raw-means.pdf", width = image.width, height = image.height, 
  pointsize = pointsize)
par(mar = mymargin)
plot(range(x[, 1]), range(x[, 2]), xlab = expression(mu[1]), ylab = expression(mu[2]), 
  main = "Raw MCMC output", type = "n")
for (k in 1:K) {
  points(mcmc.pars[thinning, k, 1:2], col = pal(K)[k], cex = 0.7, pch = 16)
}
dev.off()
# Reordered outputs for the means (9 plots)
for (method in rownames(summary(ls$permutations))) {
  pdf(file = paste("bivariate-", method, ".pdf", sep = ""), width = image.width, 
    height = image.height, pointsize = pointsize)
  par(mar = mymargin)
  plot(range(x[, 1]), range(x[, 2]), xlab = expression(mu[1]), ylab = expression(mu[2]), 
    main = method, type = "n")
  for (k in 1:K) {
    points(permute.mcmc(mcmc.pars, ls$permutations[[method]])$output[thinning, 
      k, 1:2], col = pal(K)[k], cex = 0.7, pch = 16)
  }
  dev.off()
}

# log-likelihood contour plot and clustered datapoints (9 plots)
mymargin <- c(3, 3, 2, 0.5)
pointsize <- 20
for (method in rownames(summary(ls$permutations))) {
  # estimate ergodic means
  iter <- 1
  mcmc <- permute.mcmc(mcmc.pars, ls$permutations[[method]])$output
  pars <- mcmc[iter, , ]
  for (iter in 2:iterations) {
    pars <- ((iter - 1) * pars + mcmc[iter, , ])/iter
  }
  nPoints <- 1000
  xgrid <- ygrid <- seq(-5, 5, length = nPoints)
  zValues <- array(data = NA, dim = c(nPoints, nPoints))
  for (i in 1:nPoints) {
    feta <- array(data = NA, dim = c(nPoints, 2))
    feta[, 1] <- rep(xgrid[i], nPoints)
    feta[, 2] <- ygrid
    zValues[i, ] <- normDens(feta, pars)
  }
  pdf(file = paste("bivariate-Contour-", method, ".pdf", sep = ""), width = image.width, 
    height = image.height, pointsize = pointsize)
  par(mar = mymargin)
  contour(xgrid, ygrid, zValues, lty = "solid", nlevels = 20, main = method)
  points(x, pch = 16, cex = 0.7, col = pal(K)[ls$clusters[method, ]])
  dev.off()
}

# real density and clusters
realPars <- cbind(t(mu.real), 1/t.real[1, 1, ], 1/t.real[2, 2, ], t.real[1, 2, ], 
  p.real/sum(p.real))
nPoints <- 1000
xgrid <- ygrid <- seq(-5, 5, length = nPoints)
zValues <- array(data = NA, dim = c(nPoints, nPoints))
for (i in 1:nPoints) {
  feta <- array(data = NA, dim = c(nPoints, 2))
  feta[, 1] <- rep(xgrid[i], nPoints)
  feta[, 2] <- ygrid
  zValues[i, ] <- normDens(feta, realPars)
}
pdf(file = "bivariate-Contour-True.pdf", width = image.width, height = image.height, 
  pointsize = pointsize)
par(mar = mymargin)
contour(xgrid, ygrid, zValues, lty = "solid", nlevels = 20, main = "TRUE")
points(x, pch = 16, cex = 0.7, col = pal(K)[z.real])
dev.off()

# Example with simulated dataset 2 (K = 9)
# Simulation set-up for true parameter values
# K: number of components
# n: number of observations
# d: dimensionality of multivariate normal distribution
# x: data
# mu.real: real values for the means
# t.real: real values for inverse covariance matrices
# z.real: real allocation of observations among the K mixture components
# p.real: real values for the weights of the mixture
set.seed(123123)
K <- 9
colors <- brewer.pal(K, name = "Set1")
pal <- colorRampPalette(colors)
d <- 2
x <- array(data = NA, dim = c(n, d))
n <- 280
myDistance <- 6

p.real <- rep(1, K)
p.real[K] <- 2
mu.real <- array(data = NA, dim = c(d, K))
K1 <- K - 1
tr <- seq(0, 2 * pi, length = K1 + 1)
tr <- tr[-(K1 + 1)]
mu.real[1, 1:K1] <- cos(tr)
mu.real[2, 1:K1] <- sin(tr)
mu.real[1, K] <- 0
mu.real[2, K] <- 0
mu.real <- myDistance * mu.real
t.real <- array(data = 0, dim = c(d, d, K))
for (k in 1:K) {
  for (j in 1:d) t.real[j, j, k] <- 1
}
for (j in 1:d) t.real[j, j, K] <- 1/4
if (d > 1) {
  t.real[1, 2, 1] <- 0
  t.real[2, 1, 1] <- 0
}
z.real <- sample(1:K, n, prob = p.real, replace = TRUE)
x <- array(data = NA, dim = c(n, d))
for (k in 1:K) {
  index <- which(z.real == k)
  diaspora <- solve(t.real[, , k])
  x[index, ] <- rmvnorm(length(index), mean = mu.real[, k], sigma = diaspora, method = c("eigen", 
    "svd", "chol"), pre0.9_9994 = FALSE)
}
plot(x, col = pal(K)[z.real])  #scatterplot of simulated data with true cluster labels  
# run gibbs sampler
iterations <- 20000
burn <- 5000
gs <- gibbsSampler(iterations = iterations, K = 9, x = x, burn = burn)  # [WARNING]: run-time ~2 hours   
# so better save the results when run finishes: e.g.: R> save(gs,file =
# 'example2.RData')

iterations <- iterations - burn
zChain <- gs$z
mcmc.pars <- gs$mcmc
mapindex <- gs$MLindex
allocProbs <- gs$p
gs <- 0
# apply user-defined constraint to: mu_1 - 2mu_2
newMCMC <- array(data = NA, dim = c(iterations, K, 7))
newMCMC[, , 1:6] <- mcmc.pars
for (k in 1:K) {
  newMCMC[, k, 7] <- mcmc.pars[, k, 1] - 2 * mcmc.pars[, k, 2]
}
newConstraint <- aic(newMCMC, constraint = 7)

# run label.switching.  In order to compare all methods with true clusters:
# groundTruth = z.real
set <- c("STEPHENS", "ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2", "AIC", "DATA-BASED", 
  "USER-PERM")
ls <- label.switching(method = set, zpivot = zChain[mapindex, ], z = zChain, K = K, 
  p = allocProbs, prapivot = mcmc.pars[mapindex, , ], mcmc = mcmc.pars, data = x, 
  groundTruth = z.real, userPerm = newConstraint$permutations)
# get the similarity of single best clusterings (and check the similarity of
# estimate with the ground truth and among them)
ls$similarity

# get a random sample of 500 mcmc draws for plotting
thinning <- sample(1:dim(mcmc.pars)[1], size = 500, replace = FALSE)
image.width <- 6
image.height <- 6
pointsize <- 20
mymargin <- c(4, 4, 2, 0.5)

# raw MCMC output for the means (with label switching)
pdf(file = "bivariate-raw-means-2.pdf", width = image.width, height = image.height, 
  pointsize = pointsize)
par(mar = mymargin)
plot(c(-6.5, 6.5), c(-6.5, 6.5), xlab = expression(mu[1]), ylab = expression(mu[2]), 
  main = "Raw MCMC output", type = "n")
for (k in 1:K) {
  points(mcmc.pars[thinning, k, 1:2], col = pal(K)[k], cex = 0.7, pch = 16)
}
dev.off()
# Reordered outputs for the means (7 plots)
for (method in rownames(summary(ls$permutations))) {
  pdf(file = paste("bivariate-", method, "-2.pdf", sep = ""), width = image.width, 
    height = image.height, pointsize = pointsize)
  par(mar = mymargin)
  plot(c(-6.5, 6.5), c(-6.5, 6.5), xlab = expression(mu[1]), ylab = expression(mu[2]), 
    main = method, type = "n")
  for (k in 1:K) {
    points(permute.mcmc(mcmc.pars, ls$permutations[[method]])$output[thinning, 
      k, 1:2], col = pal(K)[k], cex = 0.7, pch = 16)
  }
  dev.off()
}

mymargin <- c(3, 3, 2, 0.5)
pointsize <- 20
for (method in rownames(summary(ls$permutations))) {
  # estimate ergodic means
  iter <- 1
  mcmc <- permute.mcmc(mcmc.pars, ls$permutations[[method]])$output
  pars <- mcmc[iter, , ]
  for (iter in 2:iterations) {
    pars <- ((iter - 1) * pars + mcmc[iter, , ])/iter
  }
  
  nPoints <- 1000
  xgrid <- seq(-8, 8, length = nPoints)
  ygrid <- seq(-9, 9, length = nPoints)
  zValues <- array(data = NA, dim = c(nPoints, nPoints))
  for (i in 1:nPoints) {
    feta <- array(data = NA, dim = c(nPoints, 2))
    feta[, 1] <- rep(xgrid[i], nPoints)
    feta[, 2] <- ygrid
    zValues[i, ] <- normDens(feta, pars)
  }
  pdf(file = paste("bivariate-Contour-", method, "-2.pdf", sep = ""), width = image.width, 
    height = image.height, pointsize = pointsize)
  par(mar = mymargin)
  contour(xgrid, ygrid, zValues, method = "flattest", nlevels = 20, main = method, 
    lty = "solid")
  points(x, pch = 16, cex = 0.7, col = pal(K)[ls$clusters[method, ]])
  dev.off()
}
# real density and clusters
realPars <- cbind(t(mu.real), 1/t.real[1, 1, ], 1/t.real[2, 2, ], t.real[1, 2, ], 
  p.real/sum(p.real))
nPoints <- 1000
xgrid <- seq(-8, 8, length = nPoints)
ygrid <- seq(-9, 9, length = nPoints)
zValues <- array(data = NA, dim = c(nPoints, nPoints))
for (i in 1:nPoints) {
  feta <- array(data = NA, dim = c(nPoints, 2))
  feta[, 1] <- rep(xgrid[i], nPoints)
  feta[, 2] <- ygrid
  zValues[i, ] <- normDens(feta, realPars)
}
pdf(file = "bivariate-Contour-True-2.pdf", width = image.width, height = image.height, 
  pointsize = pointsize)
par(mar = mymargin)
contour(xgrid, ygrid, zValues, lty = "solid", nlevels = 20, main = "TRUE")
points(x, pch = 16, cex = 0.7, col = pal(K)[z.real])
dev.off()
