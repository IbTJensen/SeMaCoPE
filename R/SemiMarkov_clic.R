# SemiMarkov with CLIC-based model selection -----------------------------------
# A variant of SemiMarkov() that keeps the structure, signature and return
# value of the original but changes how a model is chosen from the
# (R_within, R_between, sat_within, sat_between) grid.
#
# WHY
# ---
# SemiMarkov() selects by the raw maximised composite log-likelihood. That is
# only a fair comparison when the candidate models carry the same optimism
# (expected in-sample minus out-of-sample fit). Across `sat` they do not: the
# nominal parameter count is identical at every saturation, but `sat` acts as
# a coarsening parameter on the interaction covariates, so smaller values give
# the same parameters a coarser design to explain and let them absorb more
# noise. On simulated data with the true interaction parameters set to zero --
# where every `sat` is equally (ir)relevant -- raw-likelihood selection chose
# sat <= 3 in 7/10 (p = 2) and 9/10 (p = 3) replicates.
#
# CLIC (Varin & Vidoni 2005; the composite-likelihood analogue of Takeuchi's
# TIC) replaces the parameter count with the effective number of parameters
#
#   edf  = tr( J H^-1 ),    H = sensitivity,  J = variability of the score
#   CLIC = -2 * cl(betahat) + 2 * edf        (smaller is better)
#
# Under a true likelihood with independent contributions J == H and edf == p,
# recovering AIC. Here the composite likelihood sums spatially correlated
# terms, so J > H and edf > p -- and, critically, edf varies with `sat` even
# though p does not. Both matrices are already computed by
# Standard_error_matrix(): H is `Sensitivity` and J is
# `Sensitivity + Second_order_term`.
#
# WHAT ELSE CHANGES
# -----------------
# Only fits that both (a) converged and (b) have a non-singular sensitivity
# matrix are eligible for selection. SemiMarkov() applies neither filter, so
# non-converged fits -- which cluster at small `sat`, where the coarse design
# invites quasi-complete separation and coefficients run to 20+ -- were being
# compared against converged ones as equals.
#
# COST
# ----
# Materially slower than SemiMarkov(). The original evaluates only the
# log-likelihood per grid point; CLIC additionally needs H and J, so
# Standard_error_matrix() runs at every grid point rather than once at the end.
#
# NOTE ON THE RETURNED `w`
# ------------------------
# The final Standard_error_matrix() call is made exactly as in SemiMarkov(),
# passing S$w directly. Standard_error_matrix() adds a `lambda` column to that
# data.table by reference, so the returned `w` carries it, as it does with
# SemiMarkov(). This is preserved deliberately for compatibility rather than
# silently changed.

# Fit one grid point and return the pieces needed for CLIC selection.
# Never throws: a fit that errors, fails to converge, or has a singular
# sensitivity matrix comes back flagged rather than aborting the whole grid.
clic_internal <- function(
  X, covariate, edgecorrection, R_within,
  R_between, sat_within, sat_between, standardize, Poisson,
  rcond_tol = 1e-10
) {
  bad <- function(reason) {
    list(
      loglik = NA_real_, edf = NA_real_, clic = NA_real_,
      convergence = NA_integer_, converged = FALSE, singular = NA,
      status = reason
    )
  }

  S <- tryCatch(
    SemiMarkov_fixed_R(
      X,
      covariate = covariate,
      edgecorrection = edgecorrection,
      R_within = R_within,
      R_between = R_between,
      sat_within = sat_within,
      sat_between = sat_between,
      standardize = standardize,
      Poisson = Poisson,
      quiet = TRUE
    ),
    error = function(e) e
  )
  if (inherits(S, "error")) {
    return(bad(paste("fit failed:", conditionMessage(S))))
  }

  loglik <- S$maximum_log_likelihood
  conv <- S$convergence
  converged <- isTRUE(conv == 0)

  # Standard_error_matrix() modifies w by reference, so hand it a copy.
  se <- tryCatch(
    Standard_error_matrix(
      X = X,
      w = data.table::copy(S$w),
      betahat = S$betahat,
      R_within = R_within,
      R_between = R_between,
      sat_within = sat_within,
      sat_between = sat_between
    ),
    error = function(e) e
  )
  if (inherits(se, "error")) {
    # solve(S) inside Standard_error_matrix() is what fails on a singular
    # sensitivity matrix, so this is the usual route to singular == TRUE.
    return(list(
      loglik = loglik, edf = NA_real_, clic = NA_real_,
      convergence = conv, converged = converged, singular = TRUE,
      status = paste("SE failed:", conditionMessage(se))
    ))
  }

  H <- se$Sensitivity
  rc <- tryCatch(rcond(H), error = function(e) NA_real_)
  if (!is.finite(rc) || rc < rcond_tol) {
    return(list(
      loglik = loglik, edf = NA_real_, clic = NA_real_,
      convergence = conv, converged = converged, singular = TRUE,
      status = sprintf("sensitivity matrix ill-conditioned (rcond = %.3g)", rc)
    ))
  }

  edf <- tryCatch(
    sum(diag((H + se$Second_order_term) %*% solve(H))),
    error = function(e) NA_real_
  )
  if (!is.finite(edf)) {
    return(list(
      loglik = loglik, edf = NA_real_, clic = NA_real_,
      convergence = conv, converged = converged, singular = TRUE,
      status = "effective df not computable"
    ))
  }

  list(
    loglik = loglik,
    edf = edf,
    clic = -2 * loglik + 2 * edf,
    convergence = conv,
    converged = converged,
    singular = FALSE,
    status = if (converged) "ok" else "did not converge"
  )
}

#' Fits semi-parametric Markov model, selecting the interaction ranges and
#' saturation parameters by CLIC rather than by raw log composite likelihood.
#'
#' Identical to [SemiMarkov()] except for how a model is chosen from the grid:
#' candidates are ranked by the composite likelihood information criterion
#' `-2 * loglik + 2 * tr(J H^-1)` (smaller is better), and only fits that
#' converged and have a non-singular sensitivity matrix are eligible.
#'
#' @inheritParams SemiMarkov
#' @param rcond_tol Reciprocal-condition-number threshold below which the
#' sensitivity matrix is treated as singular and the fit is excluded from
#' selection.
#' @return The same list as [SemiMarkov()], with the `likelihoods` element
#' extended: it now carries `loglik`, `edf`, `clic`, `convergence`,
#' `converged`, `singular` and `status` for every grid point, plus an
#' `eligible` flag and the selected row marked in `selected`.
#' @export
SemiMarkov_clic <- function(
  X, covariate, edgecorrection = NULL, R_within, R_between,
  sat_within = Inf, sat_between = Inf, standardize = TRUE,
  Poisson = FALSE, ncores = 1, quiet = FALSE, rcond_tol = 1e-10
) {
  if (
    length(R_within) == 1 &
      length(R_between) == 1 &
      length(sat_within) == 1 &
      length(sat_between) == 1
  ) {
    S <- SemiMarkov_fixed_R(
      X,
      covariate = covariate,
      edgecorrection = edgecorrection,
      R_within = R_within,
      R_between = R_between,
      sat_within = sat_within,
      sat_between = sat_between,
      standardize = standardize,
      Poisson = Poisson,
      quiet = quiet
    )

    if (!isTRUE(S$convergence == 0) && !quiet) {
      warning(
        "SemiMarkov_clic(): the single fitted model did not converge: ",
        S$convergence_message
      )
    }

    std_err <- Standard_error_matrix(
      X = X,
      w = S$w,
      betahat = S$betahat,
      R_within = R_within,
      R_between = R_between,
      sat_within = sat_within,
      sat_between = sat_between
    )

    out <- c(S, std_err)
    return(out)
  }

  opt <- expand.grid(
    R_within = R_within,
    R_between = R_between,
    sat_within = sat_within,
    sat_between = sat_between
  )
  R_omit <- opt$R_within == opt$R_between
  opt <- opt[!R_omit, ]

  eval_one <- function(i) {
    clic_internal(
      X,
      covariate = covariate,
      edgecorrection = edgecorrection,
      R_within = opt[i, "R_within"],
      R_between = opt[i, "R_between"],
      sat_within = opt[i, "sat_within"],
      sat_between = opt[i, "sat_between"],
      standardize = standardize,
      Poisson = Poisson,
      rcond_tol = rcond_tol
    )
  }

  if (ncores == 1) {
    res <- vector("list", nrow(opt))
    for (i in 1:nrow(opt)) {
      if (!quiet) {
        cat(paste("Fitting model", i, "out of", nrow(opt)), "\r")
      }
      res[[i]] <- eval_one(i)
    }
  }

  if (ncores > 1) {
    future::plan(future::multisession, workers = ncores)
    res <- furrr::future_map(
      1:nrow(opt),
      eval_one,
      .progress = TRUE
    )
  }

  liks <- data.table::data.table(
    opt,
    loglik = vapply(res, function(z) as.numeric(z$loglik), numeric(1)),
    edf = vapply(res, function(z) as.numeric(z$edf), numeric(1)),
    clic = vapply(res, function(z) as.numeric(z$clic), numeric(1)),
    convergence = vapply(res, function(z) as.integer(z$convergence), integer(1)),
    converged = vapply(res, function(z) isTRUE(z$converged), logical(1)),
    singular = vapply(res, function(z) isTRUE(z$singular), logical(1)),
    status = vapply(res, function(z) as.character(z$status), character(1))
  )
  liks[, eligible := converged & !singular & is.finite(clic)]

  if (!any(liks$eligible)) {
    stop(
      "No grid point is eligible for selection: of ", nrow(liks),
      " fits, ", sum(!liks$converged), " did not converge and ",
      sum(liks$singular, na.rm = TRUE), " had a singular sensitivity matrix. ",
      "Inspect the grid, or widen it towards larger saturation parameters."
    )
  }

  elig <- which(liks$eligible)
  w <- elig[which.min(liks$clic[elig])]
  liks[, selected := seq_len(.N) == w]

  if (!quiet) {
    n_drop <- sum(!liks$eligible)
    if (n_drop > 0) {
      message(
        "SemiMarkov_clic(): ", n_drop, " of ", nrow(liks),
        " grid points excluded from selection (",
        sum(!liks$converged), " not converged, ",
        sum(liks$singular, na.rm = TRUE), " singular)."
      )
    }
    w_ll <- elig[which.max(liks$loglik[elig])]
    if (w_ll != w) {
      message(
        "SemiMarkov_clic(): CLIC and raw log-likelihood disagree. ",
        "CLIC selects row ", w, " (sat_within = ", liks$sat_within[w],
        ", sat_between = ", liks$sat_between[w], "); raw likelihood would ",
        "select row ", w_ll, " (sat_within = ", liks$sat_within[w_ll],
        ", sat_between = ", liks$sat_between[w_ll], ")."
      )
    }
  }

  S <- SemiMarkov_fixed_R(
    X,
    covariate = covariate,
    edgecorrection = edgecorrection,
    R_within = opt[w, "R_within"],
    R_between = opt[w, "R_between"],
    sat_within = opt[w, "sat_within"],
    sat_between = opt[w, "sat_between"],
    standardize = standardize,
    Poisson = Poisson,
    quiet = quiet
  )

  std_err <- Standard_error_matrix(
    X,
    S$w,
    S$betahat,
    opt[w, "R_within"],
    opt[w, "R_between"],
    opt[w, "sat_within"],
    opt[w, "sat_between"]
  )
  out <- c(S, std_err, likelihoods = list(liks))
  return(out)
}
