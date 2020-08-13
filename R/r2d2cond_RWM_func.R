SampleBetaRWM <- function(Beta, Beta.var, n, p, X, y, sigma2, a, b, alpha.0, 
                          beta.0, log.post) {
  # Sample Beta using Random Walk Metropolis.
  # Returns Beta, whether it was accepted, and the new log posterior.
  Beta.can <- stats::rnorm(p, Beta, Beta.var)
  log.post.can <- CalcLogPost(X, y, Beta.can, sigma2, n, p, a, b, alpha.0, 
                              beta.0)
  acc.prob <- exp(log.post.can - log.post)
  accept <- (stats::runif(1) < acc.prob)
  if (accept) {
    Beta <- Beta.can
    log.post <- log.post.can
  }
  return(list(Beta = Beta, accept = accept, log.post = log.post))  
}

SampleSigma2RWM <- function(sigma2, sigma2.var, X, y, Beta, n, p, a, b, 
                            alpha.0, beta.0, log.post) {
  # Sample sigma^2 using Random Walk Metropolis.
  # Returns sigma2, whether it was accepted, and the new log posterior.
  log.sigma2.can <- stats::rnorm(1, log(sigma2), sigma2.var)
  sigma2.can <- exp(log.sigma2.can)
  log.post.can <- CalcLogPost(X, y, Beta, sigma2.can, n, p, a, b, 
                              alpha.0, beta.0)
  acc.prob <- exp(log.post.can - log.post)
  accept <- (stats::runif(1) < acc.prob)
  if (accept) {
    sigma2 <- sigma2.can 
    log.post <- log.post.can 
  }
  return(list(sigma2 = sigma2, accept = accept, log.post = log.post))  
}
