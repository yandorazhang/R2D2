#' r2d2cond
#'
#' \code{r2d2cond} adopts MCMC sampling algorithms, aiming to obtain a sequence
#' of posterior samples based on the conditional distribution on R-squared.
#'
#' @param x Scaled design matrix with no intercept (such that diag(x'x)=1).
#' @param y Centered response variable (mean 0).
#' @param a,b Shape parameters of Beta distribution for R-Squared.
#' @param a1,b1 Shape and scale parameters of Inverse Gamma distribution on
#'   sigma^2.
#' @param nu Shape parameter of Gamma distribution for local variance
#'   components.
#' @param mu Rate parameter of Gamma distribution for local variance components
#'   (optional).
#' @param q Number of expected non-zero coefficients (if nu was not specified
#'   directly).
#' @param mcmc.n Number of MCMC iterations.
#' @param gibbs Boolean, whether to use Gibbs sampling. Applicable only for
#'   uniform-on-ellipsoid prior (nu=mu=q=NULL), and a <= p/2. Will be used by
#'   default if possible.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item{beta: Matrix (mcmc.n * p) of posterior samples.}
#'   \item{sigma2: Vector (mcmc.n) of posterior samples.}
#'   \item{theta: Vector (mcmc.n) of posterior samples, if applicable.}
#'   \item{Lambda: Matrix (mcmc.n * p) of posterior samples, if applicable.}
#'   \item{w.samples: Vector (mcmc.n) of posterior samples, if Gibbs sampler
#'   was used.}
#'   \item{z.samples: Vector (mcmc.n) of posterior samples, if Gibbs sampler
#'   was used.}
#'   \item{acc.rates: List of acceptance rates for sigma2, theta and beta.
#'  }
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item{The uniform-on-ellipsoid prior will be used if nu, mu, and q are all
#'   NULL.}
#'   \item{To fit the local shrinkage model, nu or q must be specified, and mu
#'   is optional (a suitable value will be chosen based on nu, p, and n).}
#' }
#'
#' @example inst/examples/r2d2cond.R
#'
#' @export
r2d2cond <- function(x, y, a, b, a1 = .01, b1 = .01, nu = NULL, mu = NULL,
                     q = NULL, mcmc.n = 10000, gibbs = NULL) {

  cat('Fitting model with a =', a, ', b =', b, '\n')

  # Setup and validate input
  n <- nrow(x)
  p <- ncol(x)
  if (!isTRUE(all.equal(mean(y), 0, check.attributes = F)))
    stop ("y must be centered")
  ypy <- sum(y^2)
  XpX <- t(x) %*% x
  if (!isTRUE(all.equal(diag(XpX), rep(1, p), check.attributes = F)))
    stop ("X must be standardized such that diag(XpX) = (1, ..., 1)")
  ypX <- as.vector(t(x) %*% y)

  # Eigen decomposition of XpX/n
  eigXpX <- eigen(XpX / n, symmetric = T)
  V <- eigXpX$vectors
  Dvec <- eigXpX$values

  # High dimension or not
  if (p < n) {
    cat("Fitting low-dimensional model since p < n.\n")
    hd <- FALSE

    # For gibbs
    XpXinvHalf <- V %*% diag(1 / sqrt(Dvec)) / sqrt(n) # left sqrt matrix
    XpXinv <- XpXinvHalf %*% t(XpXinvHalf)

    # For sparse model
    Dhalf = diag(sqrt(Dvec))
    Dhalfinv = diag(1 / sqrt(Dvec))
    VDhalf = V %*% Dhalf
    VDhalfinv = V %*% Dhalfinv
  } else {
    cat("Fitting high-dimensional model since p >= n.\n")
    hd <- TRUE

    r1 <- n - 1   # assume X has rank n before centering
    r2 <- p - r1
    V1 <- V[, 1:r1]
    V2 <- V[, (r1 + 1):p]
    D1 <- Dvec[1:r1] # D2 = (0, ..., 0)
    V1Dhalf    <- V1 %*% diag(sqrt(D1))
    V1Dhalfinv <- V1 %*% diag(1 / sqrt(D1))
    VDhalfinv  <- V %*% diag(c(1 / sqrt(D1), rep(1, r2)))  # not really
    VDhalf     <- V %*% diag(c(sqrt(D1), rep(1, r2)))

    # For gibbs
    XpXinvHalf = V1Dhalfinv / sqrt(n)
    XpXinv = tcrossprod(XpXinvHalf)
  }

  # Get nu (shape parameter) if q is specified
  if (!is.null(q)) {
    if (!is.null(nu)) {
      stop ('Both nu and q were specified.')
    } else {
      nu <- GetNu(p, q)
    }
  }

  # Fit sparse model if nu is not null
  if (is.null(nu)) {
    sparse <- FALSE
    nu <- mu <- NULL
  } else {
    sparse <- TRUE
    if (is.null(mu)) {
      mu <- GetMu(nu, n, p)
    }
  }

  # Fit with gibbs if possible
  if (!sparse & a <= p / 2) {
    if (is.null(gibbs)) {
      gibbs <- TRUE
      cat("Using Gibbs sampling to sample from posterior.\n")
    }
  } else {
    if (is.null(gibbs)) {
      gibbs <- FALSE
    } else {
      if (gibbs) stop("Can't use gibbs here...")
    }
  }

  # Record keeping
  Beta.samples   <- matrix(NA, mcmc.n, p)
  sigma2.samples <- rep(NA, mcmc.n)
  theta.samples  <- rep(NA, mcmc.n)
  mu.samples     <- rep(NA, mcmc.n)
  Lambda.samples <- matrix(NA, mcmc.n, p)
  w.samples      <- rep(NA, mcmc.n) # used with Gibbs
  z.samples      <- rep(NA, mcmc.n) # used with Gibbs

  # Initial values (Common)
  betahat <- as.vector(XpXinv %*% ypX) # beats lm() in hd
  Beta <- ifelse(is.na(betahat), stats::rnorm(p), betahat)
  XBeta = x %*% Beta
  meanR2 = a / (a + b)
  theta = meanR2 / (1 - meanR2)
  sigma2 <- sum(XBeta^2) / (theta * n)

  # Counters to compute acceptance rates
  att <- tmpatt <- 0
  acc.beta <- tmpacc.beta <- 0
  acc.s2 <- tmpacc.s2 <- 0
  acc.th <- tmpacc.th <- 0

  ##***************************************************************************##
  ##***************************************************************************##
  if (gibbs) {
    # Initial values
    z <- w <- 1

    for (iter in 1:mcmc.n) {
      # Sample z & w if Gibbs model
      zw <- SampleZW(x, Beta, p, n, a, b, w, sigma2)
      z <- zw$z
      w <- zw$w
      z.samples[iter] <- z
      w.samples[iter] <- w

      # Sample Beta
      Beta <- SampleBetaGibbs(z, w, betahat, XpXinvHalf, sigma2)
      Beta.samples[iter,] <- Beta

      # Sample sigma2
      sigma2 <- SampleSigma2Gibbs(x, Beta, y, n, p, z, w, alpha.0 = a1, beta.0 = b1)
      sigma2.samples[iter] <- sigma2

      if (iter %% 500 == 0) {
        cat(paste('Iteration', iter, 'complete.\n'))
      }
    }
  }
  ##***************************************************************************##
  ##***************************************************************************##
  ## Adaptive Metropolis ##
  if (!gibbs && sparse) {
    # Initial values
    Lambda <- rep(nu / mu, p) # prior mean
    cat('Fitting local shrinkage model with nu =', nu, 'and initial mu =', mu,
        '\n')
    # Starting candidate variance for Metropolis random walks
    Beta.var <- .05
    sigma2.var <- 1
    theta.var <- .1
    # Adaptive metropolis for (sigma2, theta) settings
    am.start <- 500
    am.steps <- 100

    for (iter in 1:mcmc.n) {
      # Sample Beta
      att <- att + 1
      tmpatt <- tmpatt + 1
      if (hd) {
        tmp.Beta <- SampleBetaSpHD(ypX, V1Dhalfinv, theta, sigma2, Lambda, p,
                                   r1, V1Dhalf, V2, VDhalf, iter)
      } else {
        tmp.Beta <- SampleBetaSpLD(ypX, VDhalfinv, theta, sigma2, Lambda, iter)
      }
      Beta <- tmp.Beta$Beta
      Gam  <- tmp.Beta$Gam
      Beta.samples[iter,] <- Beta

      # Sample sigma2 (and theta)
      # Adaptive Metropolis
      log.post <- CalcLogPostS2Th(ypy, ypX, VDhalfinv, Gam, sigma2, theta, a,
                                  b, alpha.0 = a1, beta.0 = b1, n, p)
      s2theta <- SampleSigmaThetaAda(sigma2, theta, s2th.var, sigma2.var,
                                     log.post, iter, am.start, ypy , ypX,
                                     VDhalfinv, Gam, a, b, alpha.0 = a1, beta.0 = b1,
                                     n, p)
      accept   <- s2theta$accept
      sigma2   <- s2theta$sigma2
      theta    <- s2theta$theta
      log.post <- s2theta$log.post

      if (accept) {
        acc.s2 <- acc.s2 + 1
        tmpacc.s2 <- tmpacc.s2 + 1
      }
      sigma2.samples[iter] <- sigma2
      theta.samples[iter] <- theta

      # Update variance matrix (s2th.var) for adaptive Metropolis
      if (iter == am.start) {
        s2th.samples <- cbind(sigma2.samples, theta.samples)
        s2th.var <- stats::cov(diff(s2th.samples[1:am.start, ]))
        s2th.bar <- colMeans(s2th.samples[1:am.start, ])
      }
      if (iter > am.start & (iter %% am.steps == 0)) {
        last.n <- iter - am.steps
        s2th.samples <- cbind(sigma2.samples, theta.samples)
        s2th.var <- onlinePCA::updateCovariance(
          C = s2th.var,
          x = s2th.samples[(last.n + 1):iter, ],
          n = last.n,
          xbar = s2th.bar)
        s2th.bar <- onlinePCA::updateMean(
          xbar = s2th.bar,
          x = s2th.samples[(last.n + 1):iter, ],
          n = last.n)
      }

      # Sample Lambda via Exchange Algorithm (Murray 2006)
      Lambda <- SampleLambda(Lambda, Gam, VDhalfinv, VDhalf, V1Dhalf, V2, nu,
                             mu, hd)
      Lambda.samples[iter,] <- Lambda

      # Tuning
      if ((att %% 100 == 0) & (iter <= mcmc.n / 2)) {
        theta.var <- TuneTheta(theta.var, tmpacc.th, tmpatt)
        sigma2.var <- TuneSigma(sigma2.var, tmpacc.s2, tmpatt)
        tmpatt <- tmpacc.beta <- tmpacc.th <- tmpacc.s2 <- 0
      }

      if (iter %% 500 == 0) {
        cat(paste('Iteration', iter, 'complete.\n'))
      }
    }
  }
  ##***************************************************************************##
  ##***************************************************************************##
  # Random walk Metropolis
  if (!gibbs && !sparse){
    # Initial values
    log.post <- CalcLogPost(x, y, Beta, sigma2, n, p, a, b, alpha.0 = a1, beta.0 = b1)
    # Starting candidate variance for Metropolis random walks
    Beta.var <- .05
    sigma2.var <- 1
    theta.var <- .1

    for (iter in 1:mcmc.n) {
      # Sample Beta
      att <- att + 1
      tmpatt <- tmpatt + 1
      tmp.Beta <- SampleBetaRWM(Beta, Beta.var, n, p, x, y, sigma2, a, b,
                                alpha.0 = a1, beta.0 = b1, log.post)
      Beta <- tmp.Beta$Beta
      accept <- tmp.Beta$accept
      log.post <- tmp.Beta$log.post
      acc.beta <- acc.beta + accept
      tmpacc.beta <- tmpacc.beta + accept
      Beta.samples[iter,] <- Beta

      # Sample sigma2 (and theta)
      tmpSigma2 <- SampleSigma2RWM(sigma2, sigma2.var, x, y, Beta, n, p, a,
                                   b, alpha.0 = a1, beta.0 = b1, log.post)
      sigma2 <- tmpSigma2$sigma2
      accept <- tmpSigma2$accept
      log.post <- tmpSigma2$log.post
      if (accept) {
        acc.s2 <- acc.s2 + 1
        tmpacc.s2 <- tmpacc.s2 + 1
      }
      sigma2.samples[iter] <- sigma2

      # Tunning
      if ((att %% 100 == 0) & (iter <= mcmc.n / 2)) {
        Beta.var <-  TuneBeta(Beta.var, tmpacc.beta, tmpatt)
        sigma2.var <- TuneSigma(sigma2.var, tmpacc.s2, tmpatt)
        tmpatt <- tmpacc.beta <- tmpacc.th <- tmpacc.s2 <- 0
      }

      if (iter %% 500 == 0) {
        cat(paste('Iteration', iter, 'complete.\n'))
      }
    }
  }
  ##***************************************************************************##
  acc.rates <- list(s2 = acc.s2 / iter,
                    theta = acc.th / iter,
                    beta = acc.beta / iter)
  return (
    list(
      beta = Beta.samples,
      sigma2 = sigma2.samples,
      theta = theta.samples,
      Lambda = Lambda.samples,
      w.samples = w.samples,
      z.samples = z.samples,
      acc.rates = acc.rates
    )
  )
}
