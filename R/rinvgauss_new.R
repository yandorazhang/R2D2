#' The Inverse Gaussian Distribution
#'
#' \code{rinvgauss_new} achieves random generation from the inverse gaussian
#' distribution.
#'
#' @param n Sample size.
#' @param mean Vector of means. Default to be 1.
#' @param shape Vector of shape parameters.
#' @param dispersion Vector of dispersion parameters.
#'
#' @keywords Internal
rinvgauss_new <- function(n, mean = 1, shape = NULL, dispersion = 1) {
  if (length(n) > 1)
    n <- length(n)
  if (n < 0)
    stop("n can't be negative")
  n <- as.integer(n)
  if (n == 0)
    return(numeric(0))
  if (!is.null(shape))
    dispersion <- 1/shape
  mu <- rep_len(mean, n)
  phi <- rep_len(dispersion, n)
  r <- rep_len(0, n)
  i <- (mu > 0 & phi > 0)
  if (!all(i)) {
    r[!i] <- NA
    n <- sum(i)
  }
  phi[i] <- phi[i] * mu[i]

  repeat {
    Y <- stats::rchisq(n, df = 1)
    X1 <- 1 + phi[i]/2 * (Y - sqrt(4 * Y/phi[i] + Y^2))
    X1[which(abs(X1) < 1e-06)] = 0  # 12/05/2014
    te = 1/(1 + X1)
    te1 = unlist(lapply(te, function(x) min(x, 1)))
    # te2 = unlist(lapply(te1, function(x) max(x,0)))
    firstroot <- as.logical(stats::rbinom(n, size = 1L, prob = te1))
    if (all(!is.na(firstroot)))
      break
  }

  r[i][firstroot] <- X1[firstroot]
  r[i][!firstroot] <- 1/X1[!firstroot]
  mu * r
}
