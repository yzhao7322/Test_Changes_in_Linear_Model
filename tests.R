# Weighted-residual change tests: V (L2) and Q (maximum norm).
# Theorems 1.1--1.4, Section 3 and equations (4.2)--(4.4).
# Only cointReg is needed for automatic Andrews bandwidth selection.
# A numeric bandwidth (number of Bartlett lags) needs only base R.

.jmva_scalar <- function(x, name, lower, upper, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x)))
    stop(name, " must be ", if (integer) "an integer " else "a number ",
         "in [", lower, ", ", upper, "].")
}

# Bartlett estimator: gamma_0 + sum_l (1-l/(L+1))*(gamma_l+gamma_l').
# Every gamma_l uses denominator N, as in (3.1).
jmva_lrv <- function(scores, bandwidth = "and") {
  scores <- as.matrix(scores)
  if (!is.numeric(scores) || any(!is.finite(scores)) || nrow(scores) < 3L)
    stop("scores must be a finite numeric matrix with at least three rows.")
  N <- nrow(scores)
  if (is.character(bandwidth)) {
    if (!identical(bandwidth, "and")) stop("Use bandwidth='and' or a numeric lag count.")
    if (!requireNamespace("cointReg", quietly = TRUE))
      stop("Install cointReg: install.packages('cointReg'); or supply a numeric bandwidth.")
    bandwidth <- cointReg::getBandwidth(scores, bandwidth = "and", kernel = "ba")
  }
  .jmva_scalar(bandwidth, "bandwidth", 0, N - 1)
  L <- as.integer(round(bandwidth))
  z <- sweep(scores, 2, colMeans(scores), "-")
  D <- crossprod(z) / N
  if (L > 0L) for (ell in seq_len(L)) {
    G <- crossprod(z[seq.int(ell + 1L, N), , drop = FALSE],
                   z[seq_len(N - ell), , drop = FALSE]) / N
    D <- D + (1 - ell / (L + 1)) * (G + t(G))
  }
  list(covariance = (D + t(D)) / 2, bandwidth = L)
}

# Reusable Monte Carlo draws from the limiting distributions.
# kappa < .5: d independent Brownian bridges on [0,1].
# kappa > .5: max of TWO independent endpoint limits Z(kappa).
# Brownian time inversion: sup_{t>=1} ||W(t)||/t^kappa has the same law as
# sup_{0<s<=1} s^(kappa-1)*||W(s)||. This avoids truncating time at a finite T.
jmva_reference <- function(d, kappa = 0, statistic = c("V", "Q"),
                           reps = 1999L, grid = 2000L, seed = NULL) {
  statistic <- match.arg(statistic)
  .jmva_scalar(d, "d", 1, Inf, TRUE)
  .jmva_scalar(kappa, "kappa", 0, 1)
  .jmva_scalar(reps, "reps", 99, Inf, TRUE)
  .jmva_scalar(grid, "grid", 100, Inf, TRUE)
  if (kappa == 0.5) stop("At kappa=.5 use Vostrikova or Darling-Erdos calibration.")
  # Seeding the reference simulation does not alter the caller's RNG stream.
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  norm_rows <- if (statistic == "V") function(z) sqrt(rowSums(z^2)) else
    function(z) apply(abs(z), 1, max)
  times <- if (kappa < 0.5) seq_len(grid) / grid else
    exp(seq(log(1e-8), 0, length.out = grid))
  dt <- diff(c(0, times))
  one_path <- function() {
    W <- matrix(rnorm(grid * d), grid, d) * sqrt(dt)
    W <- vapply(seq_len(d), function(j) cumsum(W[, j]), numeric(grid))
    if (kappa < 0.5) {
      B <- W - outer(times, W[grid, ])
      # Omit t=1, where the bridge and denominator both vanish.
      idx <- seq_len(grid - 1L)
      max(norm_rows(B)[idx] / (times[idx] * (1 - times[idx]))^kappa)
    } else max(norm_rows(W) * times^(kappa - 1))
  }
  values <- numeric(reps)
  for (i in seq_len(reps)) values[i] <- if (kappa < 0.5) one_path() else
    max(one_path(), one_path())
  structure(list(values = values, d = d, kappa = kappa,
                 statistic = statistic, reps = reps, grid = grid,
                 min_time = if (kappa > 0.5) min(times) else 0),
            class = "jmva_reference")
}

# Upper-tail approximation (4.4). Intended for small significance levels.
.jmva_vos_tail <- function(x, d, N) {
  if (x <= sqrt(d)) return(1)
  T <- 2 * log(N)
  exp(d * log(x) - x^2 / 2 - (d / 2) * log(2) - lgamma(d / 2)) *
    (T - d * T / x^2 + 4 / x^2)
}

jmva_test <- function(y, x, statistic = c("V", "Q"), kappa = 0.5,
                      calibration = c("auto", "vostrikova", "darling_erdos", "simulation"),
                      alpha = 0.05, bandwidth = "and", reference = NULL,
                      reps = 1999L, grid = 2000L, seed = NULL) {
  data_name <- paste(deparse(substitute(y)), "and", deparse(substitute(x)))
  statistic <- match.arg(statistic)
  calibration <- match.arg(calibration)
  .jmva_scalar(kappa, "kappa", 0, 1)
  .jmva_scalar(alpha, "alpha", .Machine$double.eps, 0.1)
  x <- as.matrix(x)
  if (!is.numeric(x) || !is.numeric(y) || !is.null(dim(y)) ||
      nrow(x) != length(y) || ncol(x) < 1L ||
      any(!is.finite(x)) || any(!is.finite(y)))
    stop("y must be a finite numeric vector; x a matching finite numeric design matrix.")
  N <- length(y)
  d <- ncol(x)
  if (N < max(20L, d + 3L)) stop("Too few observations for this design.")
  fit <- stats::lm.fit(x, y)
  if (fit$rank < d) stop("x is rank deficient. Include an intercept only once.")
  scores <- x * as.numeric(fit$residuals)
  lrv <- jmva_lrv(scores, bandwidth)
  D <- lrv$covariance
  eig <- eigen(D, symmetric = TRUE)
  if (min(eig$values) <= max(eig$values) * 1e-10)
    stop("The long-run covariance is singular or ill-conditioned; check data/design/bandwidth.")
  # Symmetric matrix inverse square root, never elementwise sqrt(solve(D)).
  D_inv_sqrt <- eig$vectors %*% diag(1 / sqrt(eig$values), nrow = d) %*% t(eig$vectors)
  k <- seq_len(N - 1L)
  # Explicit discrete convention: k scores and weight evaluated at t=k/N.
  Z <- vapply(seq_len(d), function(j) cumsum(scores[, j]), numeric(N))[k, , drop = FALSE]
  Z <- (Z - outer(k / N, colSums(scores))) / sqrt(N)
  white <- Z %*% D_inv_sqrt
  path <- if (statistic == "V") sqrt(rowSums(white^2)) else apply(abs(white), 1, max)
  path <- path / ((k / N) * (1 - k / N))^kappa
  trim <- if (kappa > 0.5) ceiling(sqrt(N)) else 1L
  eligible <- which(k >= trim & k <= N - trim)
  peak <- eligible[which.max(path[eligible])]
  cp <- k[peak]
  raw <- path[peak]
  # t1=t2=ceil(sqrt(N))/N, not the unscaled number of observations.
  scale <- if (kappa > 0.5) (trim / N)^(kappa - 0.5) else 1
  value <- raw * scale
  if (calibration == "auto") calibration <- if (kappa == 0.5) "vostrikova" else "simulation"
  if ((kappa == 0.5) != (calibration != "simulation"))
    stop("Use Vostrikova/Darling-Erdos only at kappa=.5; simulation otherwise.")
  if (calibration == "simulation") {
    if (is.null(reference)) reference <- jmva_reference(d, kappa, statistic, reps, grid, seed)
    if (!inherits(reference, "jmva_reference") || reference$d != d ||
        reference$kappa != kappa || reference$statistic != statistic ||
        length(reference$values) < 99L || any(!is.finite(reference$values)))
      stop("reference must match d, kappa, and statistic.")
    p <- (1 + sum(reference$values >= value)) / (length(reference$values) + 1)
    critical <- unname(quantile(reference$values, 1 - alpha, type = 1))
  } else if (calibration == "darling_erdos") {
    a <- sqrt(2 * log(log(N)))
    dimension <- if (statistic == "V") d else 1L
    b <- 2 * log(log(N)) + dimension / 2 * log(log(log(N))) - lgamma(dimension / 2)
    value <- a * raw - b
    rate <- if (statistic == "V") 2 else 2 * d
    p <- -expm1(-rate * exp(-value))
    critical <- log(rate) - log(-log1p(-alpha))
  } else {
    # Work on the raw, untransformed supremum for V* and Q*.
    dimension <- if (statistic == "V") d else 1L
    tail_at <- function(q) {
      prob <- min(1, max(0, .jmva_vos_tail(q, dimension, N)))
      if (statistic == "V") prob else -expm1(d * log1p(-prob))
    }
    p <- tail_at(raw)
    lower <- sqrt(dimension)
    upper <- max(8, lower + 2)
    while (tail_at(upper) > alpha) upper <- upper * 2
    critical <- uniroot(function(q) tail_at(q) - alpha, c(lower, upper))$root
  }
  label <- paste0(statistic, if (kappa > .5) "_heavy" else "", "(kappa=", signif(kappa, 4), ")")
  structure(list(statistic = setNames(value, label), p.value = p,
                 critical.value = critical, reject = p <= alpha, alpha = alpha,
                 change_point = cp, change_fraction = cp / N,
                 raw.supremum = raw, scale = scale, kappa = kappa,
                 norm = statistic, calibration = calibration, bandwidth = lrv$bandwidth,
                 covariance = D, coefficients = fit$coefficients,
                 candidate_indices = k[eligible], path = path,
                 method = paste("Weighted-residual", statistic, "test /", calibration),
                 data.name = data_name, alternative = "at least one regression coefficient changes",
                 reference = if (calibration == "simulation") reference else NULL),
            class = c("jmva_test", "htest"))
}

print.jmva_test <- function(x, ...) {
  cat(x$method, "\n")
  cat("  statistic =", signif(unname(x$statistic), 5),
      "; p-value =", signif(x$p.value, 4), "; reject =", x$reject, "\n")
  cat("  candidate break after row", x$change_point,
      sprintf("(fraction %.3f); Bartlett lags = %d\n", x$change_fraction, x$bandwidth))
  invisible(x)
}
