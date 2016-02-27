library(mvtnorm)
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
# gs <- 0

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

# Stephens
ls <- label.switching(method = "STEPHENS", 
                      zpivot = zChain[pivot, ],
                      z = zChain,
                      K = K, 
                      p = allocProbs,
                      data = x)

# PRA
ls <- label.switching(method = "PRA", 
                      zpivot = zChain[pivot, ],
                      z = zChain,
                      K = K, 
                      prapivot = mcmc.pars[pivot, , ],
                      mcmc = mcmc.pars,
                      data = x)

# ECR
ls <- label.switching(method = "ECR", 
                      zpivot = zChain[pivot, ],
                      z = zChain,
                      K = K, 
                      data = x)

# ECR-ITERATIVE-1
ls <- label.switching(method = "ECR-ITERATIVE-1", 
                      zpivot = zChain[pivot, ],
                      z = zChain,
                      K = K, 
                      data = x)

# ECR-ITERATIVE-2
ls <- label.switching(method = "ECR-ITERATIVE-2", 
                      zpivot = zChain[pivot, ],
                      z = zChain,
                      K = K, 
                      p = allocProbs,
                      data = x)

# SJW
ls <- label.switching(method = "SJW",
                      zpivot = zChain[pivot, ], 
                      z = zChain,
                      K = K, 
                      p = allocProbs,
                      prapivot = mcmc.pars[pivot, , ],
                      complete = complete.bivariate.normal.loglikelihood, 
                      mcmc = mcmc.pars,
                      data = x,
                      sjwinit = pivot)
# AIC
ls <- label.switching(method = "AIC",
                      zpivot = zChain[pivot, ],
                      z = zChain,
                      K = K, 
                      p = allocProbs,
                      mcmc = mcmc.pars,
                      data = x)

# DATA-BASED
ls <- label.switching(method = "DATA-BASED",
                      zpivot = zChain[pivot, ],
                      z = zChain,
                      K = K, 
                      data = x)




