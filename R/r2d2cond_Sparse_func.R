SampleBetaSpHD <- function(ypX, V1Dhalfinv, theta, sigma2, Lambda, p, r1,
                           V1Dhalf, V2, VDhalf, iter) {
  # Sample Beta from the high-dimensional sparse model
  # Returns Beta, Gam and whether it failed.
  kap.mu <- as.vector(ypX %*% V1Dhalfinv) * sqrt(theta / sigma2)
  W11 <- t(V1Dhalf) %*% (V1Dhalf * Lambda)
  W21 <- t(V2) %*% (V1Dhalf * Lambda)
  W11inv <- chol2inv(chol(W11))
  Gam.list <- SampleFB(kap.mu, W11inv / 2)
  if (Gam.list$finish) {
    Gam1 <- Gam.list$x
    fail <- FALSE
  } else {
    if (iter < 1000) {
      # Use mean direction in place of Gam1 if early in MCMC
      print(paste("Using mean direction in place on iteration", iter))
      Gam1 <- as.vector(W11 %*% kap.mu)
      Gam1 <- Gam1 / sqrt(sum(Gam1^2))
      fail <- FALSE
    } else {
      fail <- TRUE
      Gam1 <- NA
    }
  }

  # sample gam2 using conditional MVN trick
  Z <- stats::rnorm(p, 0, Lambda)
  Gam.tmp <- as.vector(crossprod(VDhalf, Z))
  Gam2 <- as.vector(Gam.tmp[(r1 + 1):p] +
                      W21 %*% W11inv %*% (Gam1 - Gam.tmp[1:r1]))

  if (fail) {
    Beta <- Gam <- NA
  } else {
    Gam <- c(Gam1, Gam2)
    Beta1 <- as.vector(tcrossprod(Gam1, V1Dhalfinv))
    Beta2 <- as.vector(tcrossprod(Gam2, V2))  # safer way to calc: V2 %*% Gam2
    Beta <- (Beta1 + Beta2) * sqrt(theta * sigma2)
  }
  return(list(Beta = Beta, Gam = Gam, fail = fail))
}

SampleBetaSpLD <- function(ypX, VDhalfinv, theta, sigma2, Lambda, iter) {
  # Sample Beta from the low-dimensional sparse model.
  # Returns Beta, Gam and whether it failed.
  kap.mu <- as.vector(ypX %*% VDhalfinv) * sqrt(theta / sigma2)
  A <- t(VDhalfinv) %*% (VDhalfinv / Lambda) / 2
  Gam.list <- SampleFB(kap.mu, A)
  if (Gam.list$finish) {
    fail <- FALSE
    Gam <- Gam.list$x
    Beta <- VDhalfinv %*% Gam * sqrt(theta * sigma2)
  } else {
    if (iter < 1000) {
      # Use mean direction in place of Gam if early in MCMC
      print(paste("Using mean direction in place on iteration", iter))
      Gam <- kap.mu / sqrt(sum(kap.mu^2))
      Beta <- VDhalfinv %*% Gam * sqrt(theta * sigma2)
      fail <- FALSE
    } else {
      fail <- TRUE
      Gam <- NA
      Beta <- NA
    }
  }
  return(list(Beta = Beta, Gam = Gam, fail = fail))
}

SampleSigmaThetaAda <- function(sigma2, theta, s2th.var, sigma2.var, log.post,
                                iter, am.start, ypy , ypX, VDhalfinv, Gam, a, b,
                                alpha.0, beta.0, n, p) {
  # Sample sigma2 and theta using adaptive metropolis sampling.
  if ((iter > am.start) & (stats::runif(1) < .95)) {
    s2th.can <- MASS::mvrnorm(1, c(sigma2, theta), (2.38^2 / 2) * s2th.var)
  } else {
    # propose sigma2 and theta using Metropolis Random Walk
    s2th.can <- stats::rnorm(2, c(sigma2, theta), sigma2.var)
  }
  if (any(s2th.can <= 0)) {
    log.post.can <- -Inf
  } else {
    log.post.can <- CalcLogPostS2Th(ypy, ypX, VDhalfinv, Gam, s2th.can[1],
                                    s2th.can[2], a, b, alpha.0, beta.0, n, p)
  }
  acc.prob <- exp(log.post.can - log.post)
  accept <- stats::runif(1) < acc.prob
  if (accept) {
    sigma2 <- s2th.can[1]
    theta  <- s2th.can[2]
    log.post <- log.post.can
  }
  return (list(sigma2 = sigma2, theta = theta, accept = accept,
               log.post = log.post))
}

SampleLambda <- function(Lambda, Gam, VDhalfinv, VDhalf, V1Dhalf, V2, nu, mu,
                         hd) {
  # Sample Lambda vector via Exchange Algorithm (Murray 2006)
  eps <- 1e-7
  EPS <- 1E20
  Gam.0 <- as.vector(VDhalfinv %*% Gam)
  p <- length(Lambda)
  can.Lambda <- Lambda
  for (j in sample(p)) {
    can.Lambda[j] <- GIGrvg::rgig(1, lambda = nu, chi = Gam.0[j]^2, psi = 2*mu)
  }
  can.Lambda <- pmin(pmax(can.Lambda, eps), EPS)
  if (hd) {
    r1 <- dim(V1Dhalf)[2]
    W11 <- t(V1Dhalf) %*% (V1Dhalf * can.Lambda)
    W21 <- t(V2) %*% (V1Dhalf * can.Lambda)
    W11inv <- chol2inv(chol(W11))
    Gam1 <- as.vector(rBingham(1, W11inv / 2))

    Z <- stats::rnorm(p, 0, can.Lambda)
    Gam.tmp <- as.vector(crossprod(VDhalf, Z))
    Gam2 <- as.vector(Gam.tmp[(r1 + 1):p] +
                        W21 %*% W11inv %*% (Gam1 - Gam.tmp[1:r1]))

    Gam.star <- c(Gam1, Gam2)
    Gam.star <- VDhalfinv %*% Gam.star
  } else {
    A.can <- t(VDhalfinv) %*% (VDhalfinv / can.Lambda) / 2
    Gam.star <- as.vector(rBingham(1, A.can))
    Gam.star <- VDhalfinv %*% Gam.star
  }
  log.acc.prob <- as.vector((Gam.star^2 * (1/can.Lambda - 1/Lambda)) / 2)

  # log.acc.prob <- ifelse(can.Lambda > EPS, -Inf, log.acc.prob)
  acc.prob <- exp(log.acc.prob)
  accept <- stats::runif(p) < acc.prob
  Lambda <- ifelse(accept, can.Lambda, Lambda)
  return(Lambda)
}

TuneTheta <- function(theta.var, tmpacc.th, tmpatt) {
  if (tmpacc.th / tmpatt < .35) theta.var <- .8 *  theta.var
  if (tmpacc.th / tmpatt > .55) theta.var <- 1.2 * theta.var
  return (theta.var)
}
