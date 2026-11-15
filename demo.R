# Run all three files from the same folder:
#   Rscript demo.R
# or in R/RStudio: source("demo.R", chdir = TRUE)

local({
  # Locate companion scripts without changing the caller's working directory.
  frames <- sys.frames()
  source_files <- lapply(frames, function(fr) fr$ofile)
  source_files <- Filter(Negate(is.null), source_files)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_dir <- if (length(source_files)) {
    source_path <- tail(source_files, 1L)[[1L]]
    # source(..., chdir=TRUE) has already entered the file's directory.
    if (!file.exists(source_path) && file.exists(basename(source_path)))
      source_path <- basename(source_path)
    dirname(normalizePath(source_path, mustWork = TRUE))
  } else if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE))
  } else getwd()
  source(file.path(script_dir, "dgp.R"), local = TRUE)
  source(file.path(script_dir, "tests.R"), local = TRUE)
  if (!requireNamespace("cointReg", quietly = TRUE))
    stop("Please first run install.packages('cointReg') for automatic Andrews bandwidth.")

  # 1. The paper's main DGP: (intercept, slope) changes after row 125.
  dat <- dgp(N = 250, beta = c(1, 1), delta = 1,
             kstar = 0.5, error = "arma", seed = 2026)
  cat("True break: after row", dat$break_index, "\n\n")

  # 2. V* and Q*: kappa=.5, Vostrikova critical values/p-values.
  v <- jmva_test(dat$y, dat$x, statistic = "V")
  q <- jmva_test(dat$y, dat$x, statistic = "Q")
  print(v)
  print(q)

  # 3. Light and heavy weights. Small reference simulations keep the demo quick.
  # For final inference, increase reps and grid; cache and reuse references.
  results <- list(v, q)
  for (kap in c(0, 2 / 3)) for (norm in c("V", "Q")) {
    ref <- jmva_reference(d = ncol(dat$x), kappa = kap, statistic = norm,
                          reps = 999, grid = 2000, seed = 103)
    ans <- jmva_test(dat$y, dat$x, statistic = norm, kappa = kap, reference = ref)
    results[[length(results) + 1L]] <- ans
  }
  tab <- do.call(rbind, lapply(results, function(ans) data.frame(
    test = names(ans$statistic), calibration = ans$calibration,
    statistic = unname(ans$statistic), p_value = ans$p.value,
    reject_5pct = ans$reject, candidate_break = ans$change_point)))
  cat("\nComparison (a candidate break is meaningful only when H0 is rejected):\n")
  print(tab, row.names = FALSE, digits = 4)

  # 4. Null and GDP2 examples; replace dat with either to run the same tests.
  null_data <- dgp(N = 250, delta = 0, error = "iid", seed = 2027)
  rr2_data <- dgp2(N = 250, delta = 1, scenario = "c2", error = "garch", seed = 2028)
  cat("\nNull example (delta=0):\n")
  print(jmva_test(null_data$y, null_data$x))
  cat("\nRR2 four-coefficient example (c2, GARCH):\n")
  print(jmva_test(rr2_data$y, rr2_data$x))
})
