# Compiled (RcppArmadillo) reimplementation of beta_estimation() --------------
# Kept in its own file, alongside beta_estimation_fast() (the pure-R version
# that only applies fix #3, finite-difference gradient still), so all three
# can be benchmarked and validated against each other, e.g.
#   identical(beta_estimation(w, betastart, p, q, Poisson)$par,
#             beta_estimation_rcpp(w, betastart, p, q, Poisson)$par)
#
# This is the "go further" C++ option: instead of a hand-derived R gradient
# passed to optim(), the entire Newton-Raphson loop -- loglik, analytic
# gradient, and analytic Hessian (here: the Sensitivity/Fisher-information
# matrix S = -Hessian, the same quantity Standard_error_matrix() already
# computes) -- is fused into one compiled pass per iteration in
# cpp/beta_estimation_rcpp.cpp, avoiding both R's per-call overhead and the
# redundant exp()/rowsum() recomputation that separate loglik()/score()/
# sensitivity() R closures would each pay for independently. See that file's
# header for the numerical-stability details (log-sum-exp softmax, step
# halving) that a naive Newton implementation is missing.
#
# MISSING COVARIATE DATA
# ----------------------
# Covariates can legitimately be missing at some spatial locations (e.g. a
# covariate image that does not cover the whole window, or NAs in a supplied
# covariate data.frame). The R beta_estimation() absorbs this implicitly:
# in type_pred(), Lambda = sum(lambda) over a group containing a single NA is
# itself NA, so *every* prob in that group becomes NA, and the whole group
# drops out of `sum(log(prob[!is.na(prob)]))`. The effective behaviour is
# therefore complete-case deletion at the level of the observed point
# (group), not the individual candidate row -- dropping only the offending
# candidate row would silently change that group's softmax denominator, since
# the probabilities within a group must sum over all p candidate types.
#
# This wrapper reproduces that behaviour explicitly: any group containing a
# non-finite design entry is removed before the compiled solver runs, the
# remaining group ids are renumbered contiguously, and the number of dropped
# groups is reported (via a message, and in `groups_dropped` of the result).
# Because the non-finite pattern depends only on the design matrix and not on
# beta, this is done once up front rather than per iteration.
#
# The compiled routine lives in src/beta_estimation_rcpp.cpp and is built into
# the package DLL by R CMD INSTALL, so no compiler is needed at run time.
beta_estimation_rcpp <- function(
  w, betastart, p, q, Poisson,
  maxit = 100, tol = 1e-8, quiet = FALSE
) {
  Wmat <- as.matrix(w[, -(1:4), with = FALSE])

  grp_dt <- w[, .(xcoord, ycoord, type_obs)]
  grp_dt[, grp := .GRP, by = .(xcoord, ycoord, type_obs)]
  grp <- grp_dt$grp

  is_true <- w$type_obs == w$j

  # Complete-case deletion by group -- see "MISSING COVARIATE DATA" above.
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
    # Renumber the surviving groups contiguously. `grp` is non-decreasing
    # (w is sorted by type_obs, xcoord, ycoord, j), so unique() returns the
    # surviving ids in ascending order and match() maps them onto 1..G.
    grp <- match(grp, unique(grp))

    if (!quiet) {
      message(
        "beta_estimation_rcpp(): dropped ", n_dropped, " of ",
        n_groups_total, " observed points with missing covariate values ",
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

  fit <- beta_newton_cpp(
    Wmat = Wmat,
    grp = as.integer(grp),
    true_idx = as.integer(true_idx),
    beta = as.numeric(betastart),
    maxit = maxit,
    tol = tol
  )

  if (fit$convergence != 0 && !quiet) {
    warning(
      "beta_estimation_rcpp() did not converge: ", fit$message,
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
    groups_dropped = n_dropped
  )
}
