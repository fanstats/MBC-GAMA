#################################################################
###################################################################
#############        MtFA Factor Analyzer   #########################
###################################################################
#################################################################
## first, the pdf functions
logpdf.t <- function(X, mu, lambda, psi, v){
  n <- dim(X)[1]
  p <- length(psi)
  q <- ifelse(is.vector(lambda), 1,  dim(lambda)[2] )
  psilambda <- 1/psi * lambda 
  M <- diag(q) + crossprod(lambda, psilambda) 
  invM <- chol2inv(chol(M))             
  invsig <- diag(1/psi)- tcrossprod( tcrossprod(psilambda, invM), psilambda) 
  log.gamma <- lgamma((v+p)/2) - lgamma(v/2) 
  lconst <- -0.5 * p * (log(v) + log(pi)) 
  logdetsig <-  sum(log(psi)) - log(det(invM)) 
  logconst  <- log.gamma + lconst  - 0.5 * logdetsig
  X_center  <- sweep(X, 2, mu)
  quadratic_term <- rowSums(tcrossprod( X_center, invsig) * X_center)
  logquadratic <- -0.5 * (v+p) * log(1 + quadratic_term  / v)
  #####
  lpdf<-  logconst  + logquadratic
  result <- list(mahal= quadratic_term, lpdf = lpdf )
  return(result) 
}

  ##############################################################################
  
  ############################
  # mtfa.est.k: this uses k-means initialization. 
  ############################
  library(fad)
  mtfa.est.k <- function(X, K, q, tol = 1e-6, v = rep(30, K), maxiter = 500, nstart= 50){
    X= as.matrix(X)
    n = nrow(X)
    p = ncol(X)
    #kmeans
    kmeans_result <- kmeans(X, centers = K, nstart = 50)
    mu <- kmeans_result$centers
    sigma <- lapply(1:K, function(k) {
      cluster_data <- X[kmeans_result$cluster == k, ]
      cov(cluster_data)
    })
    sigma <- array(unlist(sigma), dim = c(p,p,K)) #to coerce sigma into an array.
    w <- as.vector(table(kmeans_result$cluster)) / n
    clusters <- kmeans_result$cluster
    n.clusters <- table(clusters)
    n.clusters <- as.vector(n.clusters)
    lambda <- array(dim = c(p, q, K) )
    psi <- matrix(nrow =K, ncol= p  )
    for (k in 1:K){ #cov2cor here might not be necessary.
      sigmaR <- suppressWarnings(cov2cor(sigma[,,k]))
      vars <- diag(sigma[,,k])
      start <- 1/ vars
      fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q, start= start, maxit = 500, tol = tol)) 
      sd = sqrt(vars)
      lambda[,,k] <-  lmd <-  fa$loadings 
      lambda[,,k]  <- sd * lmd
      psi[k,] <- ps <-  fa$uniquenesses
      psi[k,] <-  sd^2*ps
    }
    
    initials <- list(mean = mu, Sigma =sigma, weights = w , v=v, lambda =lambda, psi = psi)
    
    llh <- - Inf
    gamma <-  matrix(0, n, K)
    for (t in 1:maxiter){
      # E - Step
      logp  <- sapply(1:K, function(k) logpdf.t(X, mu[k,], lambda[,,k], psi[k,], v[k])$lpdf)
      logp  <- sweep(logp, 2, log(w), "+")
      maxvals <- apply(logp, 1, max)
      logp  <- sweep(logp, 1, maxvals, "-")
      wprob <- exp(logp)
      sums  <- rowSums(wprob)
      gamma <- sweep(wprob, 1, sums , "/")
      u <- sapply(1:K, function(k)  
        (v[k]+p) /  (v[k] + logpdf.t(X, mu[k,], lambda[,,k], psi[k,], v[k])$mahal))
      
      # CM - Step 1
      gsum <- colSums(gamma)
      w <- gsum / n
      weights <- gamma * u
      weightsum <- colSums(weights)
      mu <- crossprod(weights, X) / weightsum
      
      # CM - Step 2: compute lambdas and psis
      for (k in 1:K){
        centered_X <- sweep(X, 2, mu[k,], "-")  # X - mu[k,]
        weighted_centered_X <- weights[,k] * centered_X
        sigma[,,k] <- crossprod(weighted_centered_X, centered_X) / weightsum[k]
        
        sigmaR <- suppressWarnings(cov2cor(sigma[,,k]))
        vars <- diag(sigma[,,k])
        start <- 1 / vars
        fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q, start= start, maxit = 500, tol = tol)) 
        sd = sqrt(vars)
        lambda[,,k] <-  lam <-  fa$loadings 
        lambda[,,k] <- sd * lam 
        psi[k,] <- ps <-  fa$uniquenesses
        psi[k,] <- sd^2*ps
        
        v_func <- function(vnew) {
          vv1 <- -digamma(vnew / 2) + log(vnew / 2) + 1 
          vv2 <- (1 / gsum[k]) * sum(gamma[,k] * (log(u[,k]) - u[,k]) )
          vv3 <- digamma((v[k] + p) / 2) - log((v[k] + p) / 2)
          return(vv1 + vv2 + vv3 )
        }
        v_min = 0.0001
        v_max = 200
        f_min <- try(v_func(v_min))
        f_max <- try(v_func(v_max))
        
        if (is.finite(f_min) && is.finite(f_max)) {
          if (f_min * f_max > 0) {
            v[k] <- ifelse( abs(f_min) < abs(f_max), v_min, v_max )
          } 
          else {
            root <- tryCatch(uniroot(v_func, interval = c(v_min, v_max))$root, 
                             error = function(e) NA)
            if (is.finite(root)) {
              v[k] <- root
            } else {
              return(error_message("Root-finding failed"))
            }
          }
        } else {
          return(error_message("Invalid boundary evaluations for roots of v_func"))
        }
        
      }
      llh0 <- llh
      ####################
      logp    <- sapply(1:K, function(k) logpdf.t(X, mu[k,], lambda[,,k], psi[k,], v[k])$lpdf)
      logp    <- sweep(logp, 2, log(w),  "+")
      maxvals <- apply(logp, 1, max)
      logp    <- sweep(logp, 1, maxvals, "-")
      wprob   <- exp(logp)
      sums    <- rowSums(wprob)
      llh     <- sum(maxvals, log(sums))
      # stopping criterion. 
      if  ( abs(llh-llh0) < tol * abs(llh) ){
        break
      }
    }
    
    gamma <- sweep(wprob, 1, sums , "/")
    clusters <- apply(gamma, 1, function(x) which.max(x) )
    iter <- t
    nparam <-  (2*K - 1) + 2 * K * p + K * (p * q - q * (q - 1)/2)
    BIC <- -2*llh +  nparam * log(n) 
    
    analyzer = list(clusters = clusters,  weights = w, v= v, u=u,
                    lambda= lambda, psi=psi, sigma = sigma, niter =t,
                    means = mu, loglik= llh, BIC = BIC,  
                    nfactors = q, converged = t < maxiter, gamma = gamma, 
                    initial_param = initials
    )
    return(analyzer)
    }  
  ############################
  
  
  ############################
  # mtfa.est.r: this uses random initialization.  
  ############################
  library(fad)
  mtfa.est.r <- function(X, K, q , v  = rep(30,K), tol =1e-6, maxiter = 500, nstart= 100, init=NULL){
    X= as.matrix(X)
    n= nrow(X)
    p= ncol(X)
    if (is.null(init)) {
      set.seed(12321)
    } else {
      set.seed(init)
    }
    
    mu <- matrix(rnorm(K*p), K)
    lambda <- array(rnorm(p*q*K), dim = c(p, q, K) )
    psi <- matrix(runif(K*p, 0.3,0.9) , nrow = K, ncol = p  )  
    sigma <- array(dim = c(p, p, K) )
    w<- rep(1/K, K)
    
    initials <- list(mean = mu, Sigma =sigma, weights = w , lambda =lambda, psi = psi, v=v)
    
    llh <- - Inf
    gamma <- u <-  matrix(0, n, K)
    for (t in 1:maxiter){
      # E - Step
      logp  <- sapply(1:K, function(k) logpdf.t(X, mu[k,], lambda[,,k], psi[k,], v[k])$lpdf)
      logp  <- sweep(logp, 2, log(w), "+")
      maxvals <- apply(logp, 1, max)
      logp  <- sweep(logp, 1, maxvals, "-")
      wprob <- exp(logp)
      sums  <- rowSums(wprob)
      gamma <- sweep(wprob, 1, sums , "/")
      u <- sapply(1:K, function(k)  
              (v[k]+p) /  (v[k] + logpdf.t(X, mu[k,], lambda[,,k], psi[k,], v[k])$mahal))
      
      # CM - Step 1
      gsum <- colSums(gamma)
      w <- gsum / n
      weights <- gamma * u
      weightsum <- colSums(weights)
      mu <- crossprod(weights, X) / weightsum
      
      # CM - Step 2: compute lambdas and psis
      for (k in 1:K){
        centered_X <- sweep(X, 2, mu[k,], "-")  # X - mu[k,]
        weighted_centered_X <- weights[,k] * centered_X
        sigma[,,k] <- crossprod(weighted_centered_X, centered_X) / weightsum[k]
        
        sigmaR <- suppressWarnings(cov2cor(sigma[,,k]))
        vars <- diag(sigma[,,k])
        start <- 1 / vars
        fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q, start= start, maxit = 500, tol = tol)) 
        sd = sqrt(vars)
        lambda[,,k] <-  lam <-  fa$loadings 
        lambda[,,k] <- sd * lam 
        psi[k,] <- ps <-  fa$uniquenesses
        psi[k,] <- sd^2*ps
        
        v_func <- function(vnew) {
          vv1 <- -digamma(vnew / 2) + log(vnew / 2) + 1 
          vv2 <- (1 / gsum[k]) * sum(gamma[,k] * (log(u[,k]) - u[,k]) )
          vv3 <- digamma((v[k] + p) / 2) - log((v[k] + p) / 2)
          return(vv1 + vv2 + vv3 )
        }
        v_min = 0.0001
        v_max = 200
        f_min <- try(v_func(v_min))
        f_max <- try(v_func(v_max))
        
        if (is.finite(f_min) && is.finite(f_max)) {
          if (f_min * f_max > 0) {
            v[k] <- ifelse( abs(f_min) < abs(f_max), v_min, v_max )
          } 
          else {
            root <- tryCatch(uniroot(v_func, interval = c(v_min, v_max))$root, 
                             error = function(e) NA)
            if (is.finite(root)) {
              v[k] <- root
            } else {
              return(error_message("Root-finding failed"))
            }
          }
        } else {
          return(error_message("Invalid boundary evaluations for roots of v_func"))
        }
        
      }
      
      llh0 <- llh
      ####################
      logp    <- sapply(1:K, function(k) logpdf.t(X, mu[k,], lambda[,,k], psi[k,], v[k])$lpdf)
      logp    <- sweep(logp, 2, log(w),  "+")
      maxvals <- apply(logp, 1, max)
      logp    <- sweep(logp, 1, maxvals, "-")
      wprob   <- exp(logp)
      sums    <- rowSums(wprob)
      llh     <- sum(maxvals, log(sums))
      # stopping criterion. 
      if  ( abs(llh-llh0) < tol * abs(llh) ){
        break
      }
    }
    
    gamma <- sweep(wprob, 1, sums , "/")
    clusters <- apply(gamma, 1, function(x) which.max(x) )
    iter <- t
    
    nparam <-  (2*K - 1) + 2 * K * p + K * (p * q - q * (q - 1)/2)
    BIC <- -2*llh +  nparam * log(n) 
    
    analyzer = list(clusters = clusters,  weights = w, v= v, u=u,
                    lambda= lambda, psi=psi, sigma = sigma, niter =t,
                    means = mu, loglik= llh, BIC = BIC,  
                    nfactors = q, converged = t < maxiter, gamma = gamma, 
                    initial_param = initials
    )
    return(analyzer)
  }
  
  ######################################################################################################################
  #######################################################################################################
  # mtfa wrapper
  #########################################
  
  ##################
  # MtFA wrapper.new: mtfa with multiple initials and emEM 
  ##################
  mtfad <- function(Y, K, q , v = rep(30, K), tol = 1e-6, maxiter = 500, nstart= 20, init = NULL, innerNStart= 3){
    if (is.data.frame(Y)) Y <- as.matrix(Y)
    llhs <- numeric( (nstart) )  
    res0 <- try({mtfa.est.k(Y, K, q,  v = rep(30, K), tol =tol,  maxiter = 1000) }, TRUE)
    llhs0 <- tryCatch({res0$loglik} , error= function(e) -Inf ) 
    
    if (is.null(init)) {
      set.seed(1234321)
    } else {
      set.seed(init)
    }  
    seed <- runif(nstart)*10^8
    for (j in 1:nstart){
      res <- try({ mtfa.est.r(Y, K, q,  v = rep(30, K), tol =tol, init= seed[j], maxiter = 5) }, TRUE)
      llhs[j] <- tryCatch({res$loglik},  error= function(e) -Inf )     
    }
    llhs1 <- c(llhs, llhs0)
    if (sum(llhs1==-Inf)== (nstart + 1)) stop("All initial values and kmeans encountered errors. 
                      Number of factors or number of clusters might have 
                        been chosen too high compared to sample size
                        Also, try standardizing data. \n")
    irank <- order(llhs, decreasing = TRUE)[1:min(innerNStart, max(1,floor(nstart/4)))]
    
    for (j in irank){
      res2 <- try({
        mtfa.est.r(Y, K, q,  v = rep(30, K),  tol =tol, init= seed[j], maxiter = maxiter) 
      }, TRUE)
      llhs[j] <- tryCatch({res2$loglik},  error= function(e) -Inf)     
    }
    llhs1[1:nstart] <- llhs
    bestindex<- which.max(llhs1)
    if (bestindex== (nstart + 1) ){
      return( res0 )
    } else {
      res.f <-  mtfa.est.r(Y, K, q,  v = rep(30, K), tol = tol, init = seed[bestindex],  maxiter = maxiter)
      res.f$init = seed[bestindex]
      return( res.f )
    }
  }
  
################################################################################
  
  
  
  
  #########################################################################################################
  #######################################################################################################
  #############################################################################
 
  
  
  
  #############################################################
  ###################################################################
  #############        Generalized MtFA Factor Analyzer : MtFAD-q   ####################
  ###################################################################
  #############################################################
  
  ############################
  # Generalized mtfa.est.kq: this uses kmeans initialization. 
  ############################
  mtfa.est.kq <- function(X, K, q, tol = 1e-6, v = rep(30, K), maxiter = 500, nstart= 50){
    X= as.matrix(X)
    n = nrow(X)
    p = ncol(X)
    #kmeans
    kmeans_result <- kmeans(X, centers = K, nstart = 50)
    mu <- kmeans_result$centers
    sigma <- lapply(1:K, function(k) {
      cluster_data <- X[kmeans_result$cluster == k, ]
      cov(cluster_data)
    })
    sigma <- array(unlist(sigma), dim = c(p,p,K)) #to coerce sigma into an array.
    w <- as.vector(table(kmeans_result$cluster)) / n
    clusters <- kmeans_result$cluster
    n.clusters <- table(clusters)
    n.clusters <- as.vector(n.clusters)
    lambda <- list() # array(dim = c(p, q, K) )
    psi <- matrix(nrow =K, ncol= p  )
    for (k in 1:K){ 
      sigmaR <- suppressWarnings(cov2cor(sigma[,,k]))
      vars <- diag(sigma[,,k])
      start <- 1/ vars
      fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q[k], start= start, maxit = 500, tol = tol)) 
      sd = sqrt(vars)
      lambda[[k]] <-  lmd <-  fa$loadings 
      lambda[[k]]  <- sd * lmd
      psi[k,] <- ps <-  fa$uniquenesses
      psi[k,] <-  sd^2*ps
    }
    
    initials <- list(mean = mu, Sigma =sigma, weights = w , v=v, lambda =lambda, psi = psi)
    
    llh <- - Inf
    gamma <-  matrix(0, n, K)
    for (t in 1:maxiter){
      # E - Step
      logp  <- sapply(1:K, function(k) logpdf.t(X, mu[k,], lambda[[k]], psi[k,], v[k])$lpdf)
      logp  <- sweep(logp, 2, log(w), "+")
      maxvals <- apply(logp, 1, max)
      logp  <- sweep(logp, 1, maxvals, "-")
      wprob <- exp(logp)
      sums  <- rowSums(wprob)
      gamma <- sweep(wprob, 1, sums , "/")
      u <- sapply(1:K, function(k)  
        (v[k]+p) /  (v[k] + logpdf.t(X, mu[k,], lambda[[k]], psi[k,], v[k])$mahal))
      # CM - Step 1
      gsum <- colSums(gamma)
      w <- gsum / n
      weights <- gamma * u
      weightsum <- colSums(weights)
      mu <- crossprod(weights, X) / weightsum
      
      # CM - Step 2: compute lambdas and psis
      for (k in 1:K){
        centered_X <- sweep(X, 2, mu[k,], "-")  # X - mu[k,]
        weighted_centered_X <- weights[,k] * centered_X
        sigma[,,k] <- crossprod(weighted_centered_X, centered_X) / weightsum[k]
        
        sigmaR <- suppressWarnings(cov2cor(sigma[,,k]))
        vars <- diag(sigma[,,k])
        start <- 1 / vars
        fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q[k], start= start, maxit = 500, tol = tol)) 
        sd = sqrt(vars)
        lambda[[k]] <-  lam <-  fa$loadings 
        lambda[[k]] <- sd * lam 
        psi[k,] <- ps <-  fa$uniquenesses
        psi[k,] <- sd^2*ps
        
        v_func <- function(vnew) {
          vv1 <- -digamma(vnew / 2) + log(vnew / 2) + 1 
          vv2 <- (1 / gsum[k]) * sum(gamma[,k] * (log(u[,k]) - u[,k]) )
          vv3 <- digamma((v[k] + p) / 2) - log((v[k] + p) / 2)
          return(vv1 + vv2 + vv3 )
        }
        v_min = 0.0001
        v_max = 200
        f_min <- try(v_func(v_min))
        f_max <- try(v_func(v_max))
        
        if (is.finite(f_min) && is.finite(f_max)) {
          if (f_min * f_max > 0) {
            v[k] <- ifelse( abs(f_min) < abs(f_max), v_min, v_max )
          } 
          else {
            root <- tryCatch(uniroot(v_func, interval = c(v_min, v_max))$root, 
                             error = function(e) NA)
            if (is.finite(root)) {
              v[k] <- root
            } else {
              return(error_message("Root-finding failed"))
            }
          }
        } else {
          return(error_message("Invalid boundary evaluations for roots of v_func"))
        }
      }
      llh0 <- llh
      ####################
      logp    <- sapply(1:K, function(k) logpdf.t(X, mu[k,], lambda[[k]], psi[k,], v[k])$lpdf)
      logp    <- sweep(logp, 2, log(w),  "+")
      maxvals <- apply(logp, 1, max)
      logp    <- sweep(logp, 1, maxvals, "-")
      wprob   <- exp(logp)
      sums    <- rowSums(wprob)
      llh     <- sum(maxvals, log(sums))
      ## stopping criterion. 
      if  ( abs(llh-llh0) < tol * abs(llh) ){
        break
      }
    }
    
    gamma <- sweep(wprob, 1, sums , "/")
    clusters <- apply(gamma, 1, function(x) which.max(x) )
    iter <- t
    nparam <-  nparam <- (2*K - 1) + 2 * K * p + sum(p * q - q * (q - 1)/2)
    BIC <- -2*llh +  nparam * log(n) 
    
    analyzer = list(clusters = clusters,  weights = w, v= v, u=u,
                    lambda= lambda, psi=psi, sigma = sigma, niter =t,
                    means = mu, loglik= llh, BIC = BIC,  
                    nfactors = q, converged = t < maxiter, gamma = gamma, 
                    initial_param = initials
    )
    return(analyzer)
  }  
  ############################
  
  
  ############################
  #Generalized  mtfa.est.rq: this uses random initialization.  
  ############################
  mtfa.est.rq <- function(X, K, q , v  = rep(30,K), tol =1e-6, maxiter = 500, nstart= 100, init=NULL){
    X= as.matrix(X)
    n= nrow(X)
    p= ncol(X)
    if (is.null(init)) {
      set.seed(12321)
    } else {
      set.seed(init)
    }
    
    mu <- matrix(rnorm(K*p), K)
    lambda <- list() #array(rnorm(p*q*K), dim = c(p, q, K) )
    psi <- matrix(runif(K*p, 0.3,0.9) , nrow = K, ncol = p  )  
    sigma <- array(dim = c(p, p, K) )
    w<- rep(1/K, K)
    
    initials <- list(mean = mu, Sigma =sigma, weights = w , lambda =lambda, psi = psi, v=v)
    
    llh <- - Inf
    gamma <- u <-  matrix(0, n, K)
    for (t in 1:maxiter){
      # E - Step
      logp  <- sapply(1:K, function(k) logpdf.t(X, mu[k,], lambda[[k]], psi[k,], v[k])$lpdf)
      logp  <- sweep(logp, 2, log(w), "+")
      maxvals <- apply(logp, 1, max)
      logp  <- sweep(logp, 1, maxvals, "-")
      wprob <- exp(logp)
      sums  <- rowSums(wprob)
      gamma <- sweep(wprob, 1, sums , "/")
      u <- sapply(1:K, function(k)  
        (v[k]+p) /  (v[k] + logpdf.t(X, mu[k,], lambda[[k]], psi[k,], v[k])$mahal))
      
      # CM - Step 1
      gsum <- colSums(gamma)
      w <- gsum / n
      weights <- gamma * u
      weightsum <- colSums(weights)
      mu <- crossprod(weights, X) / weightsum
      
      # CM - Step 2: compute lambdas and psis, and v's
      for (k in 1:K){
        centered_X <- sweep(X, 2, mu[k,], "-")  
        weighted_centered_X <- weights[,k] * centered_X
        sigma[,,k] <- crossprod(weighted_centered_X, centered_X) / weightsum[k]
        
        sigmaR <- suppressWarnings(cov2cor(sigma[,,k]))
        vars <- diag(sigma[,,k])
        start <- 1 / vars
        fa <- suppressWarnings(fad:::fad.fit.cor(R = sigmaR, q= q[k], start= start, maxit = 500, tol = tol)) 
        sd = sqrt(vars)
        lambda[[k]] <-  lam <-  fa$loadings 
        lambda[[k]] <- sd * lam 
        psi[k,] <- ps <-  fa$uniquenesses
        psi[k,] <- sd^2*ps
        
        v_func <- function(vnew) {
          vv1 <- -digamma(vnew / 2) + log(vnew / 2) + 1 
          vv2 <- (1 / gsum[k]) * sum(gamma[,k] * (log(u[,k]) - u[,k]) )
          vv3 <- digamma((v[k] + p) / 2) - log((v[k] + p) / 2)
          return(vv1 + vv2 + vv3 )
        }
        v_min = 0.0001
        v_max = 200
        f_min <- try(v_func(v_min))
        f_max <- try(v_func(v_max))
        
        if (is.finite(f_min) && is.finite(f_max)) {
          if (f_min * f_max > 0) {
            v[k] <- ifelse( abs(f_min) < abs(f_max), v_min, v_max )
          } 
          else {
            root <- tryCatch(uniroot(v_func, interval = c(v_min, v_max))$root, 
                             error = function(e) NA)
            if (is.finite(root)) {
              v[k] <- root
            } else {
              return(error_message("Root-finding failed"))
            }
          }
        } else {
          return(error_message("Invalid boundary evaluations for roots of v_func"))
        }
        
      }
      
      llh0 <- llh
      ####################
      logp    <- sapply(1:K, function(k) logpdf.t(X, mu[k,], lambda[[k]], psi[k,], v[k])$lpdf)
      logp    <- sweep(logp, 2, log(w),  "+")
      maxvals <- apply(logp, 1, max)
      logp    <- sweep(logp, 1, maxvals, "-")
      wprob   <- exp(logp)
      sums    <- rowSums(wprob)
      llh     <- sum(maxvals, log(sums))
      # stopping criterion. 
      if  ( abs(llh-llh0) < tol * abs(llh) ){
        break
      }
    }
    
    gamma <- sweep(wprob, 1, sums , "/")
    clusters <- apply(gamma, 1, function(x) which.max(x) )
    iter <- t
    
    nparam <-  (2*K - 1) + 2 * K * p + sum(p * q - q * (q - 1)/2)
    BIC <- -2*llh +  nparam * log(n) 
    
    analyzer = list(clusters = clusters,  weights = w, v= v, u=u,
                    lambda= lambda, psi=psi, sigma = sigma, niter =t,
                    means = mu, loglik= llh, BIC = BIC,  
                    nfactors = q, converged = t < maxiter, gamma = gamma, 
                    initial_param = initials
    )
    return(analyzer)
  }
  ##################
  # Generalized MtFA wrapper.new: mtfad.q with multiple initials and emEM 
  ##################
  mtfad.q <- function(Y, K, q , v = rep(30, K), tol = 1e-6, maxiter = 500, nstart= 20, init = NULL, innerNStart= 3){
    if (is.data.frame(Y)) Y <- as.matrix(Y)
    if (length(q)!= K ) stop("Length of q is not compatible with K. \n")
    llhs <- numeric( (nstart) )  
    res0 <- try({mtfa.est.kq(Y, K, q, v = rep(30, K), tol =tol,  maxiter = 1000) }, TRUE)
    llhs0 <- tryCatch({res0$loglik} , error= function(e) -Inf ) 
    
    if (is.null(init)) {
      set.seed(1234321)
    } else {
      set.seed(init)
    }  
    seed <- runif(nstart)*10^8
    for (j in 1:nstart){
      res <- try({ mtfa.est.rq(Y, K, q, v = rep(30, K),tol =tol, init= seed[j], maxiter = 5) }, TRUE)
      llhs[j] <- tryCatch({res$loglik},  error= function(e) -Inf )     
    }
    llhs1 <- c(llhs, llhs0)
    if (sum(llhs1==-Inf)== (nstart + 1)) stop("All initial values and kmeans encountered errors. 
                      Number of factors or number of clusters might have 
                        been chosen too high compared to sample size
                        Also, try standardizing data. \n")
    irank <- order(llhs, decreasing = TRUE)[1:min(innerNStart, max(1,floor(nstart/4)))]
    
    for (j in irank){
      res2 <- try({
        mtfa.est.rq(Y, K, q, v = rep(30, K), tol =tol, init= seed[j], maxiter = maxiter) 
      }, TRUE)
      llhs[j] <- tryCatch({res2$loglik},  error= function(e) -Inf)     
    }
    llhs1[1:nstart] <- llhs
    bestindex<- which.max(llhs1)
    if (bestindex== (nstart + 1) ){
      return( res0 )
    } else {
      res.f <-  mtfa.est.rq(Y, K, q, v = rep(30, K), tol = tol, init = seed[bestindex],  maxiter = maxiter)
      res.f$init = seed[bestindex]
      return( res.f )
    }
  }
  ################################################################################
  
  ################################################################################
  
  ################################################################################
  
  ################################################################################
  get_init_para<- function(p,g,q,init=NULL, v= rep(30, g) ){
    K=g
    if (is.null(init)) {
      set.seed(12321)
    } else {
      set.seed(init)
    }
    
    mu <- matrix(rnorm(K*p),K )
    lambda <- array(rnorm(p*q*K), dim = c(p, q, K) )
    psi <- matrix(runif(K*p, 0.3,0.9) , nrow =K, ncol= p  )  
    w<- rep(1/K, K)
    psiarray <- array(dim=c(p,p, K))
    for (k in 1:K){psiarray[,,k]<- diag(psi[k,])}
    init_para <- list( g = K, q = q, pivec = w, mu  = t(mu), B = lambda, D = psiarray,  
                       sigma_type = "unique", D_type = "unique", v=v)
    return(init_para)
  }
  
  emmix_t_emEM <- function(Y, K, q, tol = 1e-6, maxiter = 500, nstart= 20, df_init= rep(30, K), init = NULL, innerNStart= 3){
    if (is.data.frame(Y)) Y <- as.matrix(Y)
    llhs <- numeric( (nstart) )  

    
    if (is.null(init)) {
      set.seed(123456)
    } else {
      set.seed(init)
    }  
    seed <- runif(nstart)*10^8
    for (j in 1:nstart){
      init_para = get_init_para(p, K, q, v= df_init, init=seed[j])
      res <- mtfa(Y, g=K, q=q, sigma_type = "unique", D_type = "unique", tol = tol, 
                 conv_measure = 'ratio', itmax = 5, init_para= init_para)
      llhs[j] <- tryCatch({res$logL},  error= function(e) -Inf )     
    }
 
    if (sum(llhs==-Inf)== (nstart)) stop("All initial values encountered errors. 
                      Number of factors or number of clusters might have 
                        been chosen too high compared to sample size
                        Also, try standardizing data. \n")
    irank <- order(llhs, decreasing = TRUE)[1:min(innerNStart, max(1,floor(nstart/4)))]
    for (j in irank){
      init_para = get_init_para(p,K,q,init=seed[j])
      res2 <- mtfa(Y, g=K, q=q, sigma_type = "unique", D_type = "unique", tol = tol, 
                  conv_measure = 'ratio', itmax = maxiter, init_para= init_para)
      llhs[j] <- tryCatch({res2$logL},  error= function(e) -Inf )      
    }
    bestindex<- which.max(llhs)
    init_para = get_init_para(p,K,q,init=seed[bestindex])
    res.f <- mtfa(Y, g=K, q=q, sigma_type = "unique", D_type = "unique", tol = tol, 
                 conv_measure = 'ratio', itmax = maxiter, init_para= init_para)
    
    return( res.f )
  }
  #############################
  
  
  
  
  
  
  
  
  