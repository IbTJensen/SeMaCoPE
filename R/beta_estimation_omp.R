# Multi-threaded (OpenMP) reimplementation of beta_estimation() ----------------
# Parallel counterpart of beta_estimation_rcpp(); kept in its own file so the
# serial and threaded versions can be benchmarked and validated against each
# other, e.g.
#   all.equal(beta_estimation_rcpp(w, betastart, p, q, Poisson)$par,
#             beta_estimation_omp (w, betastart, p, q, Poisson)$par)
#
# Same signature as beta_estimation()/beta_estimation_rcpp(), plus `nthreads`.
# The missing-data handling below is identical to beta_estimation_rcpp() --
# see that file for the full rationale. In short: a covariate missing at some
# spatial location makes that observed point's Lambda (and hence every prob in
# its group) NA in the R beta_estimation(), so the whole point drops out of
# `sum(log(prob[!is.na(prob)]))`. That is complete-case deletion at the level
# of the observed point, and is reproduced explicitly here before the compiled
# solver runs.
#
# CHOOSING nthreads
# -----------------
# The parallel gain is concentrated in the O(n k^2) Sensitivity matrix, so
# threading pays off when the number of parameters k is large and is close to
# pure overhead when it is small (k <~ 50). The default follows
# getOption("SeMaCoPE.threads"), falling back to the OpenMP maximum.
#
# If R is linked against a *threaded* BLAS (OpenBLAS, MKL) rather than the
# reference BLAS, the inner dgemm/syrk calls are themselves threaded, and
# running both layers at once oversubscribes the machine. In that case either
# set nthreads = 1 here and let the BLAS parallelise, or pin the BLAS to one
# thread (e.g. RhpcBLASctl::blas_set_num_threads(1)) and let this layer do it.
# Check with sessionInfo()$BLAS.
beta_estimation_omp <- function(
  w, betastart, p, q, Poisson,
  nthreads = getOption("SeMaCoPE.threads", NULL),
  maxit = 100, tol = 1e-8, quiet = FALSE
) {
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("beta_estimation_omp() requires the Rcpp package.")
  }
  if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
    stop("beta_estimation_omp() requires the RcppArmadillo package.")
  }
  if (!exists("beta_newton_omp_cpp", mode = "function", envir = .GlobalEnv)) {
    cpp_path <- file.path("cpp", "beta_estimation_omp.cpp")
    if (!file.exists(cpp_path)) {
      stop(
        "Could not find ", cpp_path, ". beta_estimation_omp() expects to be ",
        "run with the package root as the working directory."
      )
    }
    Rcpp::sourceCpp(cpp_path)
  }

  if (is.null(nthreads)) nthreads <- omp_max_threads_cpp()
  nthreads <- max(1L, as.integer(nthreads))

  Wmat <- as.matrix(w[, -(1:4), with = FALSE])

  grp_dt <- w[, .(xcoord, ycoord, type_obs)]
  grp_dt[, grp := .GRP, by = .(xcoord, ycoord, type_obs)]
  grp <- grp_dt$grp

  is_true <- w$type_obs == w$j

  # Complete-case deletion by group (see header).
  bad_row <- rowSums(!is.finite(Wmat)) > 0
  n_groups_total <- max(grp)
  if (any(bad_row)) {
    bad_grp <- unique(grp[bad_row])
    keep <- !(grp %in% bad_grp)
    n_dropped <- length(bad_grp)

    if (!any(keep)) {
      stop(
        "All ", n_groups_total, " observed points have at least one missing ",
        "(non-finite) covariate value, leaving no data to fit. Check for a ",
        "covariate that is entirely NA, or for a zero standard deviation ",
        "column being standardized (sd == 0 gives 0/0)."
      )
    }

    Wmat <- Wmat[keep, , drop = FALSE]
    grp <- grp[keep]
    is_true <- is_true[keep]
    grp <- match(grp, unique(grp))

    if (!quiet) {
      message(
        "beta_estimation_omp(): dropped ", n_dropped, " of ", n_groups_total,
        " observed points with missing covariate values ",
        "(complete-case by point, matching beta_estimation())."
      )
    }
  } else {
    n_dropped <- 0L
  }

  true_idx <- which(is_true)
  n_groups_used <- max(grp)
  if (length(true_idx) != n_groups_used) {
    stop(
      "Expected exactly one row with type_obs == j per observed point, but ",
      "found ", length(true_idx), " such rows for ", n_groups_used,
      " points. This usually means duplicated (xcoord, ycoord, type_obs) ",
      "coordinates have merged distinct points into one group."
    )
  }

  fit <- beta_newton_omp_cpp(
    Wmat = Wmat,
    grp = as.integer(grp),
    true_idx = as.integer(true_idx),
    beta = as.numeric(betastart),
    nthreads = nthreads,
    maxit = maxit,
    tol = tol
  )

  if (fit$convergence != 0 && !quiet) {
    warning(
      "beta_estimation_omp() did not converge: ", fit$message,
      " Returned estimates are the last accepted iterate."
    )
  }

  par <- as.numeric(fit$par)
  names(par) <- names(betastart)

  list(
    par = par,
    value = fit$value,
    convergence = fit$convergence,
    message = fit$message,
    iterations = fit$iterations,
    gradient = fit$gradient,
    Sensitivity = fit$Sensitivity,
    groups_used = n_groups_used,
    groups_dropped = n_dropped,
    threads_used = fit$threads_used
  )
}
