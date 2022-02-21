## Use this file as a source for "bayesian.robit.bothalgo.R" and
## "bayesian.robit.bothalgo.high_dim.R" files.                   

#------------------------------------------------------------------------------------
## Functions for finding square root and inverse  ##
## of a positive definite matrix below:           ##
#------------------------------------------------------------------------------------

# Function to find the square root of a P.D. matrix using Cholesky decomposition.
sqrt.mat <- function(B) {
  t(chol(B))
}

# Function to find inverse of a P.D. matrix using Cholesky decomposition.
solve.chol <- function(B){
  chol2inv(chol(B))
}

#------------------------------------------------------------------------------------
## Functions for Probit regression below: ##
#------------------------------------------------------------------------------------

# Function to generate from Truncated Normal distribution
# truncated on (-Inf, 0) if Y = 0, and (0, Inf) if Y = 1.
rTN <- function(mu, sigma, omega) {
  n <- length(mu)
  omega_1 <- which(omega == 1)
  if(length(omega_1) == 0){
    return(rtruncnorm(n, a = rep(-Inf,n), b = 0, mean = mu, sd = rep(1,n)))
  } else if(length(omega_1) == n) {
    return(rtruncnorm(n, a = 0, b = rep(Inf,n), mean = mu, sd = rep(1,n)))
  } else {
    lower.all <- rep(-Inf, n)
    lower.all[omega_1] <- 0
    upper.all <- rep(Inf, n)
    upper.all[-omega_1] <- 0
    rtruncnorm(n, a = lower.all, b = upper.all, mean = mu, sd = rep(1,n))
  }
}

# The proper probit DA chain
run_probit_beta_MCMC <- function(X, y, beta.start = rnorm(ncol(X)),
                                 beta_a = rep(0, ncol(X)),
                                 Sigma_a = diag(1, ncol(X)),
                                 nsim = 1000, burn.in = 100, ...){
  n <- nrow(X)
  p <- ncol(X)
  beta <- matrix(0, nrow = nsim+1, ncol = p)
  beta[1,] <- beta.start
  beta.var <- solve(t(X)%*%X + Sigma_a)
  beta.var.v <- as.numeric(beta.var %*% (Sigma_a %*% beta_a))
  beta.var.Xt <- beta.var %*% t(X)
  
  beta.var.sqrt <- sqrt.mat(beta.var)
  beta.var.sqrt.Xt <- beta.var.sqrt %*% t(X)
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  
  for(k in 1:nsim) {
    z <- rTN(mu = as.numeric(X%*%beta[k,]), sigma = rep(1,n), omega = y)
    beta.mean <- as.numeric(beta.var.Xt %*% z)+ beta.var.v
    
    beta[(k+1),] <- beta.mean + as.numeric(beta.var.sqrt %*% rnorm(p))
    
    setTxtProgressBar(pb, value = k)
  }
  return(beta[-(1:burn.in), , drop = F])
}


# The function for drawing observations from a density
# proportional to w*(g) as described in Roy.V (2012) paper.
# The following function draws a 'single' observation from w*(g).
rwg <- function(A, B, N, eps = 0.1){
  if(A <= 0) stop("non-positive A")
  
  if(B == 0) u <- rgamma(n=1, shape=(N/2), scale=(2/A))
  else {
    M <- exp(B^2/(2*eps*A)) * gamma(N/2) * (2/((1-eps)*A))^(N/2) 
    # NOTE: For non-zero B / non-zero prior mean, we need to take  
    # eps = 0.5 in M to be consistent with u and rho defined below.
    
    u <- rgamma(n=1, shape=(N/2), scale=(4/A))
    l <- u^(N/2 - 1) * exp(-A*u/2 + B*sqrt(u))
    rho <- l/(M * dgamma(u, shape=(N/2), scale=(4/A))) 
    if(runif(1) > rho) rwg(A, B, N , eps)
  }
  sqrt(u)
}



# The proper probit sandwich chain
run_probit_sandwich_beta_MCMC <- function(X, y, beta.start = rnorm(ncol(X)),
                                          beta_a = rep(0, ncol(X)),
                                          Sigma_a = diag(1, ncol(X)),
                                          nsim = 1000, burn.in = 100, ...){
  n <- nrow(X)
  p <- ncol(X)
  
  beta.s <- matrix(0, nrow = nsim+1, ncol = p)
  beta.s[1,] <- beta.start
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  Sigma_a.beta_a <- Sigma_a %*% beta_a
  
  beta.s.var <- solve(t(X)%*%X + Sigma_a)
  beta.s.var.v <- as.numeric(beta.s.var %*% Sigma_a.beta_a)
  beta.s.var.Xt <- beta.s.var %*% t(X)
  B.mat <- X %*% beta.s.var
  B.mat.v <- B.mat %*% Sigma_a.beta_a
  A.mat <- diag(1, n, n) - B.mat %*% t(X)
  
  beta.s.var.sqrt <- sqrt.mat(beta.s.var)
  beta.s.var.sqrt.Xt <- beta.s.var.sqrt %*% t(X)
  
  v.norm.sq <- sum((as.vector(Sigma_a.beta_a))^2 )

  if(v.norm.sq == 0) {
    for(k in 1:nsim) {
      z.s <- rTN(mu=as.numeric(X%*%beta.s[k,]), sigma=rep(1,n), omega=y)
      
      A <- sum(z.s * (A.mat %*% z.s))
      B <- 0
      
      g.s <- rwg(A, B, N = n)
      z1.s <- g.s * z.s
      
      beta.s.mean <- as.numeric(beta.s.var.Xt %*% z1.s)
      beta.s[(k+1),] <- beta.s.mean + as.numeric(beta.s.var.sqrt %*% rnorm(p))
      
      setTxtProgressBar(pb, value = k)
    }
  } else {
    for(k in 1:nsim) {
      z.s <- rTN(mu=as.numeric(X%*%beta.s[k,]), sigma=rep(1,n), omega=y)
      
      A <- sum(z.s * (A.mat %*% z.s))
      B <- sum(z.s * B.mat.v)
      
      g.s <- rwg(A, B, N = n)
      z1.s <- g.s * z.s
      
      beta.s.mean <- as.numeric(beta.s.var.Xt %*% z1.s) + beta.s.var.v
      beta.s[(k+1),] <- beta.s.mean + as.numeric(beta.s.var.sqrt %*% rnorm(p))
      
      setTxtProgressBar(pb, value = k)
    }
  }
  return(beta.s[-(1:burn.in), , drop = F])
}


#------------------------------------------------------------------------------------
## Functions for Robit Regression below: ##
#------------------------------------------------------------------------------------

# Function for generating from Truncated Student's t distribution with nu
# degrees of freedom.
rTt <- function(mu, nu, y){
  n <- length(mu)
  omega_1 <- which(y == 1)
  if(length(omega_1) == 0){
    return(rtt(n, location = mu, scale = 1, df = nu, left = -Inf, right = 0))
  } else if(length(omega_1) == n) {
    return(rtt(n, location = mu, scale = 1, df = nu, left = 0, right = Inf))
  } else {
    lower.all <- rep(-Inf, n)
    lower.all[omega_1] <- 0
    upper.all <- rep(Inf, n)
    upper.all[-omega_1] <- 0
    rtt(n, location = mu, scale = 1, df = nu,
        left = lower.all, right = upper.all)
  }
}


# The proper robit DA chain
run_robit_beta_MCMC <- function(X, y, nu = 3, beta.start = rnorm(ncol(X)),
                                beta_a = rep(0, ncol(X)),
                                Sigma_a = diag(1, ncol(X)),
                                nsim = 1000, burn.in = 100, ...){
  n <- nrow(X)
  p <- ncol(X)
  beta <- matrix(0, nrow = nsim+1, ncol = p)
  beta[1,] <- beta.start
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  Sigma_a.beta_a <- Sigma_a %*% beta_a
  
  for(k in 1:nsim) {
    Xbeta <- X%*%beta[k,]
    z <- rTt(mu = as.numeric(Xbeta), nu = nu, y = y)
    lambda <- rgamma(n = n, shape = 0.5*(nu+1),
                     rate = 0.5 * (nu + (z - Xbeta)^2) )
    
    tx.lambda <- t(X) %*% diag(lambda, n)
    beta.var <- solve(tx.lambda %*% X + Sigma_a)
    
    beta.var.sqrt <- tryCatch(sqrt.mat(beta.var), error = function(x){x})
    if(is(beta.var.sqrt, "error"))browser()
    
    beta.mean <- as.numeric( beta.var %*% (tx.lambda%*% z + Sigma_a.beta_a) )
    
    beta[(k+1),] <- beta.mean + as.numeric(beta.var.sqrt %*% rnorm(p))
  
    setTxtProgressBar(pb, value = k)
  }
  return(beta[-(1:burn.in), , drop = F])
}


# The proper robit sandwich chain
run_robit_sandwich_beta_MCMC <- function(X, y, nu = 3,
                                         beta.start = rnorm(ncol(X)),
                                         beta_a = rep(0, ncol(X)),
                                         Sigma_a = diag(1, ncol(X)),
                                         nsim = 1000, burn.in = 100, ...){
  n <- nrow(X)
  p <- ncol(X)
  
  beta.s <- matrix(0, nrow = nsim+1, ncol = p)
  beta.s[1,] <- beta.start
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  Sigma_a.beta_a <- Sigma_a %*% beta_a
  
  for(k in 1:nsim) {
    Xbeta <- X %*% beta.s[k, ]
    z.s <- rTt(mu = as.numeric(Xbeta), nu = nu, y = y)
    lambda.s <- rgamma(n = n, shape = 0.5*(nu+1),
                       rate = 0.5 * (nu + (z.s - Xbeta)^2) )
    
    tx.lambda <- t(X) %*% diag(lambda.s, n)
    beta.s.var <- solve(tx.lambda %*% X + Sigma_a)
    
    # middle step
    tx.lambda.half <- t(X) %*% diag(sqrt(lambda.s), n)
    Q <- t(tx.lambda.half) %*% beta.s.var %*% tx.lambda.half
    lambda.half.z <- diag(sqrt(lambda.s), n) %*% z.s
    g.s.squared.gamma.scale <- t(lambda.half.z) %*% (diag(1, n) - Q) %*% lambda.half.z
    g.s <- sqrt(rgamma(n = 1, shape = 0.5*n,
                       rate = 0.5*as.numeric(g.s.squared.gamma.scale)) )
    z1.s <- g.s * z.s
    
    # last step
    beta.s.mean <- as.numeric( beta.s.var %*% (tx.lambda %*% z1.s + Sigma_a.beta_a) )
    beta.s.var.sqrt <- tryCatch(sqrt.mat(beta.s.var), error = function(x){x})
    beta.s[(k+1),] <- beta.s.mean + as.numeric(beta.s.var.sqrt %*% rnorm(p))
    
    setTxtProgressBar(pb, value = k)
  }

  return(beta.s[-(1:burn.in), , drop = F])
}


#------------------------------------------------------------------------------------
## Log-likelihood and log-posterior calculations       ##
## for both probit and robit regression models below:  ##
#------------------------------------------------------------------------------------

# Function to calculate Log-likelihood and log-posterior for probit regression.
log.lik.post.probit <- function(X, y, beta, beta_a = rep(0, ncol(X)),
                                Sigma_a = diag(1, ncol(X)), ...){
  n <- nrow(X)
 
  quantile.x.beta.vec <- X %*% beta
  log.Phi.x.beta.vec <- pnorm(q = quantile.x.beta.vec, log.p = T)
  log.Phi.x.beta.vec.upper <- pnorm(q = quantile.x.beta.vec, log.p = T,
                                    lower.tail = F)
  out.loglik <- sum( (y * log.Phi.x.beta.vec) + ((1-y) * log.Phi.x.beta.vec.upper) )
  
  out.log.prior <- dmvnorm(x = beta, mean = beta_a, sigma = solve(Sigma_a),
                           log = T, checkSymmetry = F)
  out.logposterior <- out.loglik + out.log.prior
  
  out <- list(logLikelihood = out.loglik, logPosterior = out.logposterior)
  return(out)
}


# Function to calculate Log-likelihood and log-posterior for robit regression.
log.lik.post.robit <- function(X, y, nu = 3, beta, beta_a = rep(0, ncol(X)),
                               Sigma_a = diag(1, ncol(X)), ...){
  n <- nrow(X)
  
  quantile.x.beta.vec <- X %*% beta
  log.F_nu.x.beta.vec <- pt(q = quantile.x.beta.vec, df = nu, log.p = T)
  log.F_nu.x.beta.vec.upper <- pt(q = quantile.x.beta.vec, df = nu, log.p = T,
                                  lower.tail = F)
  out.loglik <- sum( (y * log.F_nu.x.beta.vec) + ((1-y) * log.F_nu.x.beta.vec.upper) )
  
  out.log.prior <- dmvnorm(x = beta, mean = beta_a, sigma = solve(Sigma_a),
                           log = T, checkSymmetry = F)
  out.logposterior <- out.loglik + out.log.prior
  
  out <- list(logLikelihood = out.loglik, logPosterior = out.logposterior)
  return(out) 
}


