# DGP
# Testing for changes in linear models using weighted residuals.
# This file only defines functions. No packages are required.

# y_i = x_i' beta_i + error_i. A break k means beta changes AFTER row k.
# Scalar delta changes every coefficient; a vector selects coefficients.
dgp <- function(N = 250, beta = c(1, 1), delta = 0,
                kstar = 0.5, error = c("iid", "arma", "garch"),
                x_sd = rep(1, length(beta) - 1L),
                standardize = TRUE, burnin = 1000L, seed = NULL) {
  error <- match.arg(error)
  if (length(N) != 1L || !is.finite(N) || N != floor(N) || N < 10)
    stop("N must be an integer >= 10.")
  if (!is.numeric(beta) || !length(beta) || any(!is.finite(beta)))
    stop("beta must be a finite numeric vector, including the intercept.")
  d <- length(beta)
  if (!is.numeric(delta) || !(length(delta) %in% c(1L, d)) ||
      any(!is.finite(delta))) stop("delta must be a scalar or length(beta) vector.")
  if (length(kstar) != 1L || !is.finite(kstar) || kstar <= 0 || kstar >= 1)
    stop("kstar must be strictly between zero and one.")
  k <- floor(N * kstar)
  if (k < 1L || k >= N) stop("kstar must leave observations on both sides.")
  if (length(x_sd) != d - 1L || any(!is.finite(x_sd)) || any(x_sd <= 0))
    stop("x_sd must contain one positive standard deviation per slope.")
  if (length(burnin) != 1L || !is.finite(burnin) ||
      burnin < 0 || burnin != floor(burnin)) stop("burnin must be a nonnegative integer.")
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize))
    stop("standardize must be TRUE or FALSE.")
  if (!is.null(seed)) set.seed(seed)
  delta <- rep_len(delta, d)
  x <- matrix(1, N, d)
  if (d > 1L) x[, -1L] <- sweep(matrix(rnorm(N * (d - 1L)), N), 2, x_sd, "*")
  colnames(x) <- c("intercept", if (d > 1L) paste0("x", seq_len(d - 1L)))

  e <- switch(error,
    iid = rnorm(N),
    arma = as.numeric(stats::arima.sim(list(ar = 0.3, ma = 0.2), n = N)),
    garch = {
      # Standard GARCH(1,1), matching fGarch's model in the original scripts:
      # h_i = .01 + .05 * e_{i-1}^2 + .9 * h_{i-1}; e_i = sqrt(h_i)*eta_i.
      z <- rnorm(N + burnin)
      values <- numeric(length(z))
      h_previous <- 0.01 / (1 - 0.05 - 0.9)
      e_previous <- sqrt(h_previous) * rnorm(1)
      for (i in seq_along(z)) {
        h <- 0.01 + 0.05 * e_previous^2 + 0.9 * h_previous
        values[i] <- sqrt(h) * z[i]
        h_previous <- h
        e_previous <- values[i]
      }
      tail(values, N)
    })
  # Retain the sample-SD normalization used in the existing simulation code.
  if (standardize) e <- e / sd(e)
  beta_path <- matrix(rep(beta, each = N), N, d)
  after <- seq.int(k + 1L, N)
  beta_path[after, ] <- sweep(beta_path[after, , drop = FALSE], 2, delta, "+")
  list(x = x, y = rowSums(x * beta_path) + e, errors = e,
       beta_before = beta, beta_after = beta + delta,
       break_index = if (any(delta != 0)) k else NA_integer_,
       kstar = kstar, error = error, standardize = standardize)
}

# RR2 c1: only coefficient 2 changes; c2: coefficients 2 and 3 change.
# Regressor standard deviations are (1, .5, 2), besides the intercept.
dgp2 <- function(N = 250, delta = 0, scenario = c("c1", "c2"),
                    kstar = 0.5, error = c("iid", "arma", "garch"),
                    beta = rep(1, 4), seed = NULL, ...) {
  scenario <- match.arg(scenario)
  error <- match.arg(error)
  if (length(beta) != 4L) stop("DGP2 requires four coefficients.")
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta))
    stop("RR2 delta must be a finite scalar.")
  shift <- c(0, delta, if (scenario == "c2") delta else 0, 0)
  dgp(N = N, beta = beta, delta = shift, kstar = kstar, error = error,
      x_sd = c(1, 0.5, 2), seed = seed, ...)
}
