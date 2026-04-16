gdtrans <- function(mapped_data, mu, sigma, prob){
  # generalized distributional transformation for composite cluster Cm with K simple groups
  # mapped_data: data matrix for Cm 
  # mu, sigma: mean - p by K, covariance matrix - p by p by K for Cm
  # prob: grouping probabilities for Cm
  
  eta = prob/sum(prob)
  
  transformed_data <- array(dim = dim(mapped_data))
  n <- dim(mapped_data)[1]
  
  # transform data according to the marginal cdfs
  for(j in 1:ncol(mapped_data)){
    Fhat <- rep(0,n)
    for (k in 1:length(eta)) {
      Fhat <- Fhat + eta[k]*pnorm(mapped_data[,j], mean = mu[j,k], sd = sqrt(diag(sigma[,,k])[j]))
    }
    transformed_data[,j] <- Fhat 
  }
  return(transformed_data)
}

# G-transformation function with marginal cdfs
Gtrans.sync <- function(data, mu, sigma, prob, cl = NULL, VariableSelection = FALSE, p_threshold = 0.05, ...){
  if (!all(sapply(data, is.numeric))){
    stop("All columns/attributes of the data need to be numeric")
  }
  data_c <- gdtrans(data, mu, sigma, prob)
  data_normal <- data.frame(qnorm(data_c, mean = 0, sd = 1))
  colnames(data_normal) <- colnames(data)
  
  if(VariableSelection){
    if(is.null(class)){
      stop("class information required for variable selection!")
    }
    mat.lm.pval <- apply(data_normal, 2, function(x, gr)(anova(lm(x~gr))[5][[1]][1]), gr = as.factor(cl))
    data_normal <- data_normal[,p.adjust(mat.lm.pval, method = "BH") < p_threshold]
  }
  return(data_normal)
}
