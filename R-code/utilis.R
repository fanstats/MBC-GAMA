# generate data from GMM using R::MixSim with three given omegas
genX <- compiler::cmpfun(function(n,p,q,K,omega = c(0.001,0.01,0.05),error = F){
  prob = abs(rnorm(K)); prob = prob/sum(prob)
  if(!error){
    init = random.init(p,q,K)
    U = tau = NULL
  }else{
    init = random.init.error(p,q,K)
    U = init$U
    tau = init$tau
  }
  M = init$M
  L = init$L
  D = init$D
  
  S = array(NA,c(p,p,K))
  for (k in 1:K) {
    S[,,k] = tcrossprod(L[,,k]) + diag(U[,k]) + diag(tau) # data with measurement errors
  }
  
  # generate X1 for given overlap rates: omega / search for sigma
  rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = S)
    
  s1 = s2 = s3 <- 1
  jump <- 0.01
  
  if(rate < omega[1]){
    
      while (rate < omega[1]) {
        s1 <- s1 + jump
        rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s1*S)
      }
      r1 = rate
      
      s2 = s1
      while (rate < omega[2]) {
        s2 <- s2 + jump
        rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s2*S)
      }
      r2 = rate
      
      s3 = s2
      while (rate < omega[3]) {
        s3 <- s3 + jump
        rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s3*S)
      }
      r3 = rate
      
  }else if(rate < omega[2]){
    
    while (rate > omega[1]) {
      s1 <- s1 - jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s1*S)
    }
    r1 = rate
    
    while (rate < omega[2]) {
      s2 <- s2 + jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s2*S)
    }
    r2 = rate
    
    s3 = s2
    while (rate < omega[3]) {
      s3 <- s3 + jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s3*S)
    }
    r3 = rate
    
  }else if(rate < omega[3]){
    
    while (rate > omega[2]) {
      s2 <- s2 - jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s2*S)
    }
    r2 = rate
    
    s1 = s2
    while (rate > omega[1]) {
      s1 <- s1 - jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s1*S)
    }
    r1 = rate
    
    while (rate < omega[3]) {
      s3 <- s3 + jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s3*S)
    }
    r3 = rate
    
  }else{
    
    while (rate > omega[3]) {
      s3 <- s3 - jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s3*S)
    }
    r3 = rate
    
    s2 = s3
    while (rate > omega[2]) {
      s2 <- s2 - jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s2*S)
    }
    r2 = rate
    
    r1  = r2
    while (rate > omega[1]) {
      s1 <- s1 - jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s1*S)
    }
    r1 = rate
    
  }
  glap = c(r1,r2,r3)
  sigma = c(s1,s2,s3)
    
  A1 <- MixSim::simdataset(n = n, Pi = prob, Mu = t(M), S = s1*S)
  A2 <- MixSim::simdataset(n = n, Pi = prob, Mu = t(M), S = s2*S)
  A3 <- MixSim::simdataset(n = n, Pi = prob, Mu = t(M), S = s3*S)
  
  #X2 <- fad:::.postmdiag(matrix(rnorm(n*p),n,p), sqrt(tau)) # data with error variances
  
  return(list(X1 = list(A1$X,A2$X,A3$X), id = list(A1$id,A2$id,A3$id), #X2 = X2, 
              M = M, L = L, D = D, U = U, tau = tau, prob = prob, 
              glap = glap, sigma = sigma))
})

# generate measurement error data Y/tau from Ga(a,tau)
# p-dim vector: a: shape; tau: scale
genY <- compiler::cmpfun(function(n,p,a,tau){
  Y <- matrix(NA,n,p)
  for (j in 1:p) {
    Y[,j] <- rgamma(n, shape = a[j], scale = tau[j])
  }
  return(Y)
})

# generate data from GMM using R::MixSim with a given omega
genX.w <- compiler::cmpfun(function(n,p,q,K,omega = 0.001,
                                    init = NULL, m = NULL, error = F){
  prob = abs(rnorm(K)); prob = prob/sum(prob)
  
  if(is.null(init)){
    if(!error){
      init = random.init(p,q,K)
      M = init$M
      L = init$L
      U = init$D
      tau = NULL
    }else{
      init = random.init.error(p,q,K)
      M = init$M
      L = init$L
      U = init$U
      tau = init$tau
    }
  }else{
    M = init$M
    L = init$L
    if(!error){
      U = init$D
      tau = NULL
    }else{
      U = init$U
      tau = init$tau
    }
  }
  if(!is.null(m)){
    init$M = m
    M = m
  }
  
  if(q>0){
    S = array(NA,c(p,p,K))
    for (k in 1:K) {
      if(is.null(tau)){
        S[,,k] = tcrossprod(L[,,k]) + diag(U[,k])
      }else{
        S[,,k] = tcrossprod(L[,,k]) + diag(U[,k]) + diag(tau) 
      }
    }
  }else{
    S = init$S
  }
  
  # generate X1 for the given overlap rate-omega / search for sigma
  rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = S)
  
  s1 <- 1
  jump <- 0.001
  
  if(rate < omega){
    while (rate < omega){
      s1 <- s1 + jump
      rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s1*S)
      }
    }else{
      while (rate > omega){
        s1 <- s1 - jump
        rate <- MixSim::overlapGOM(Pi = prob, Mu = t(M), S = s1*S)
      }
    }
  
  A1 <- MixSim::simdataset(n = n, Pi = prob, Mu = t(M), S = s1*S)
  
  return(list(X1 = A1$X, id = A1$id, init = init, prob = prob,
              glap = rate, sigma = s1))
})

#logit transformation for probability
lgp <- function(prob){
  lp = rep(NA,length(prob))
  for (i in 1:length(prob)) {
    prob[i] = max(1e-10,prob[i])
    prob[i] = min(1-1e-10,prob[i])
    lp[i] = log(prob[i]/(1-prob[i]))
  }
  return(lp);
}

#transform back to prob
plg <- function(lp){
  prob = rep(NA,length(lp))
  for (i in 1:length(lp)){
    prob[i] = 1/(1+1/exp(lp[i]))
  }
  return(prob);
}

##
## Display the map of pairwise overlap measures of Maitra and Melnykov
## (JCGS, 2012)
##
## overlap.mat = matrix of total pairwise overlaps (necessarily symmetric)
## map.col = colormap for the mapping
## linescol = color for the lines drawing the squares
## map.range = range of the overlap map (default: minimum and maximum of lower
##               triangle of the matrix)
## lab.col = color of the labels (same as nrow(matrix) if provided)
## lab.cex = character size of the label
## map.cex = character size of the overlap values laid on the map
## legend.cex = character size of the legend text (does not work always)
##
## provides map of overlap values for each group of mixture model
##
## written Ranjan Maitra, Ames, IA 50011-1210, June 28, 2009
##
## modified Ranjan Maitra, Ames, IA 50011-1090, December 23, 2016.
## last modified Ranjan Maitra, Ames, IA 50011-1090, October 31, 2020.
##
## modification to bring in specifications for label color, labels, maps and
## legend character size; subsequently change default color of map and lines
##

#overlap.mat = overlap_mat

#linescol = "#1D91C0"; map.range = NULL; lab = NULL; lab.col = 1; 
#lab.cex = 5; map.cex = 5; legend.cex = 1; font = 1; scale = NULL; 
#scale.pos = -2; legend.width = 1; utf.chars = FALSE
#lab.col = sim.palette; 
#map.cex = 2.75; lab.cex = 4;
#legend.cex = 1.5; font = 2; scale.pos = 0.35; legend.width = 0.1

overlap.map <-
  function(overlap.mat, map.col = RColorBrewer::brewer.pal(name = "PuRd", n = 9), 
           linescol = "#1D91C0", map.range = NULL, lab = NULL, lab.col = 1, 
           lab.cex = 5, map.cex = 5, legend.cex = 1, font = 1, scale = NULL, 
           scale.pos = -2, legend.width = 1, utf.chars = FALSE)
  {
    oxmat <- overlap.mat
    oxmat[lower.tri(oxmat)] <- NA
    diag(oxmat) <- NA
    p <- ncol(oxmat)
    if(is.null(lab)){
      lab <- 1:p
    }
    newox <- oxmat[-p, -1]
    newox <- cbind(rbind(NA, newox), NA)[, p:1]
    
    layout(matrix(c(rep(1, 4*p^2), rep(2, 2*p), rep(3, 2*p), rep(4, 2*p), rep(4, 2*p)), nrow = 2*p, ncol = 2*p + 4 ))
    ##  layout(matrix(c(rep(1, p^2), rep(2, p), rep(3, p), rep(4, p)), nrow = p, ncol = p + 3))
    
    
    par(mar = c(0.1,0.1,0.75,0.1))
    if (is.null(map.range)) map.range <- range(newox, na.rm = T)  
    image(x = 1:p, y = 1:p, z = newox, axes = F, xlab = "", ylab = "", #col = brewer.pal(9, "GnBu"))
          col = map.col, zlim = map.range)
    text(y = 2:p, x = rep(1, p-1), labels = lab[p:2], cex = lab.cex, font = font, 
         col = lab.col[p:2])
    text(x = 2:p, y = rep(1, (p-1)), labels = lab[1:(p-1)], cex = lab.cex, font = font, 
         col = lab.col[1:(p-1)])
    
    if (!is.null(scale)){
      text(x = p-scale.pos, y = p+0.2, labels = scale, cex = map.cex)
    }
    
    
    
    for(i in 1:p) {
      for(j in i:p) {
        lines(x = c(i+0.5, i+0.5), y = c(p-j+1,p-j)+1.5, col = linescol, lwd = 0.5)
        lines(y = c(i+0.5, i+0.5), x = c(p-j+1,p-j)+1.5, col = linescol, lwd = 0.5)
      }
    }
    for(i in 2:p) {
      text(x=1:p, y = i, labels = round(newox[,i], 3),
           col = ifelse(newox[,i] < median(map.range), "black", "white"), 
           cex = map.cex)
    }
    
    frame()
    
    # savepar <- par(cex=0.75, lwd=0.25, mar = c(1, 0.5, 1, 2),
    #                 xaxs="i", yaxs="i")
    
    savepar <- par(cex=0.75, lwd=0.25, mar = c(0.5, 0.5, 0.5, 0.5),
                   xaxs="i", yaxs="i")
    # if (legend){
    plot.new()
    length.col <- length(map.col) + 1
    ra <- seq(from = map.range[1], to = map.range[2], length=length.col)
    plot.window(xlim=c(0,0.1), ylim= c(map.range[1], map.range[2]))
    rect(0, ra[-length.col], legend.width, ra[-1], col = map.col, border = NULL)
    axis(4, at = ra, labels = round(ra, digits = 2), las = 1, cex.axis = legend.cex, line = NA)
    
    
    rect(0, 1, legend.width, ra[length.col], col = NULL)
    frame()
  }

