SampleZW <- function(X, Beta, p, n, a, b, w, sigma2) {
  # Sample z and w for Gibbs model.
  BXXB <- sum((X %*% Beta)^2)
  # sample z
  tmp.alpha <- p / 2 + b
  tmp.beta  <- BXXB / (2 * w * sigma2) + n / 2
  z <- 1 / stats::rgamma(1, tmp.alpha, tmp.beta)

  # sample w
  tmp.alpha <- p / 2 - a
  tmp.beta <- BXXB / (2 * z * sigma2)
  w.theta <- stats::rgamma(1, tmp.alpha, tmp.beta)
  w <- 1 / (1 + w.theta)
  if (w < .Machine$double.eps) w <- .Machine$double.eps # a too small

  return(list(z = z, w = w))
}

SampleBetaGibbs <- function(z, w, betahat, XpXinvHalf, sigma2) {
  # Sample Beta for Gibbs model.
  C <- (z * w) / (1 + z * w) # shrinkage factor
  M <- C * betahat
  VV <- XpXinvHalf * sqrt(C * sigma2)
  return(rMvNormLeft(mu = M, SigmaHalf = VV))
}

SampleSigma2Gibbs <- function(X, Beta, y, n, p, z, w, alpha.0, beta.0) {
  # Sample sigma^2 for Gibbs model.
  XBeta <- X %*% Beta
  BXXB  <- sum((XBeta)^2)
  sse   <- sum((y - XBeta)^2)
  alpha.star <- (n + p) / 2 + alpha.0
  beta.star  <- sse / 2 + BXXB / (2 * z * w) + beta.0
  sigma2 <- 1 / stats::rgamma(1, alpha.star, beta.star)
  return(sigma2)
}
