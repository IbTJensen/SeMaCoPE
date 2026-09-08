# Covariate_setup with separate within/between saturation ----------------------
# A variant of Covariate_setup_rcpp() that replaces the single `sat` argument
# with `sat_within` and `sat_between`, applied exactly the way `R_within` and
# `R_between` already are:
#
#   R  [m,l] = R_within   if m == l, else R_between      (unchanged)
#   sat[m,l] = sat_within if m == l, else sat_between    (new)
#
# so a within-type interaction saturates at sat_within and a between-type one
# at sat_between. Either may be Inf independently, so sat_within = Inf with a
# finite sat_between gives Strauss behaviour within types and Geyer saturation
# between them.
#
# Everything else -- the covariate handling, the design-matrix layout, the
# column names and ordering, the returned list(w, q) -- is identical to
# Covariate_setup_rcpp(). Setting sat_within == sat_between == sat reproduces
# Covariate_setup_rcpp(..., sat = sat) exactly.
#
# NOT WIRED IN. SemiMarkov_fixed_R() still calls Covariate_setup_rcpp() with a
# single `sat`; nothing else in the package has been changed to accommodate
# this variant. To use it, call it directly, or thread sat_within/sat_between
# through SemiMarkov_fixed_R(), SemiMarkov() and Standard_error_matrix()
# (the latter derives its Int_range from `sat == Inf`, which would need to
# become a per-pair test).
#
# The compiled core lives in src/covariate_setup_sat2_rcpp.cpp.
#
# The mark-ordering caveat from Covariate_setup_rcpp() applies here too: this
# indexes points directly and so is correct for any input ordering, whereas
# the original R Covariate_setup() requires X to be sorted by mark.
Covariate_setup_sat2_rcpp <- function(
  X, Xis, nis, covariate, R_within,
  R_between, sat_within, sat_between, Poisson, mark.pp
) {
  p <- length(Xis)
  n <- X$n

  prelim_dt <- data.table::data.table(
    xcoord = X$x,
    ycoord = X$y,
    type_obs = as.integer(X$marks)
  )

  if (!Poisson) {
    # Map marks onto contiguous 1..p codes. For factor marks this is the
    # identity on as.integer(marks); doing it explicitly also makes integer
    # marks with arbitrary values (e.g. 2, 5, 7) work.
    type_code <- match(as.integer(X$marks), as.integer(mark.pp))
    if (anyNA(type_code)) {
      stop("Covariate_setup_sat2_rcpp(): every mark of X must appear in mark.pp.")
    }

    core <- covariate_setup_core_sat2_cpp(
      x = as.numeric(X$x),
      y = as.numeric(X$y),
      type = as.integer(type_code),
      p = p,
      R_within = R_within,
      R_between = R_between,
      sat_within = as.numeric(sat_within),
      sat_between = as.numeric(sat_between),
      want_components = FALSE
    )
    Int_mat <- core$Int_mat
    colnames(Int_mat) <- paste(rep(1:p, each = p), rep(1:p, p), sep = "-")
  }

  # Construct interaction matrix
  pts <- lapply(1:p, function(j) {
    data.table::data.table(prelim_dt[, 1:3], j = j)
  })
  pts <- data.table::rbindlist(pts)

  # Spatial covariates
  ## Intercept
  Intercept_matrix <- model.matrix(~ factor(pts$j, levels = c(p, 1:(p - 1))))
  Intercept_matrix <- Intercept_matrix[, -1, drop = F]
  colnames(Intercept_matrix) <- paste("Intercept", 1:(p - 1), sep = ":")
  spat_cov <- cbind(pts, Intercept_matrix)
  q <- 1
  ## Covariate input
  covar.function <- function(X, covariate, covar_name, spat_cov) {
    Z <- data.table::data.table(
      xcoord = X$x,
      ycoord = X$y,
      type_obs = as.integer(X$marks),
      matrix(covariate(X$x, X$y), X$n, p - 1)
    )
    colnames(Z)[4:(p + 2)] <- paste(covar_name, 1:(p - 1), sep = ":")

    spat_cov <- merge(spat_cov, Z, by = c("xcoord", "ycoord", "type_obs"))
    Intercept_cols <- colnames(spat_cov)[grepl("Intercept", colnames(spat_cov))]
    Covariate_cols <- colnames(spat_cov)[grepl(covar_name, colnames(spat_cov))]
    cov_dt <- spat_cov[, ..Covariate_cols] * spat_cov[, ..Intercept_cols]
    spat_cov[, (ncol(spat_cov) - p + 2):(ncol(spat_cov)) := cov_dt]
    return(spat_cov)
  }

  covar.im <- function(X, covariate, covar_name, spat_cov) {
    Z <- data.table::data.table(
      xcoord = X$x,
      ycoord = X$y,
      type_obs = as.integer(X$marks),
      matrix(covariate[X], X$n, p - 1)
    )
    colnames(Z)[4:(p + 2)] <- paste(covar_name, 1:(p - 1), sep = ":")

    spat_cov <- merge(spat_cov, Z, by = c("xcoord", "ycoord", "type_obs"))
    Intercept_cols <- colnames(spat_cov)[grepl("Intercept", colnames(spat_cov))]
    Covariate_cols <- colnames(spat_cov)[grepl(covar_name, colnames(spat_cov))]
    cov_dt <- spat_cov[, ..Covariate_cols] * spat_cov[, ..Intercept_cols]
    spat_cov[, (ncol(spat_cov) - p + 2):(ncol(spat_cov)) := cov_dt]
    return(spat_cov)
  }

  if (!is.null(covariate)) {
    if (spatstat.geom::is.im(covariate)) {
      spat_cov <- covar.im(X, covariate, covar_name = "Covariate", spat_cov)
      q <- 2
    }
    if (is.function(covariate)) {
      spat_cov <- covar.function(
        X,
        covariate,
        covar_name = "Covariate",
        spat_cov
      )
      q <- 2
    }
    if (is.data.frame(covariate) | is.matrix(covariate)) {
      if (ncol(covariate) < 3) {
        stop(
          "When spatial covariates are given as a matrix/data.frame, ",
          "it must contain columns with the x- and y-coordinates and at ",
          "least one covariate."
        )
      }
      if (is.null(colnames(covariate))) {
        stop(
          "When spatial covariates are given as a matrix/data.frame, ",
          "the first two columns must contain the x- and y-coordinates, ",
          "with \nnames 'xcoord' and 'ycoord'."
        )
      }
      if (!all(colnames(covariate)[1:2] == c("xcoord", "ycoord"))) {
        stop(
          "When spatial covariates are given as a matrix/data.frame, ",
          "the first two columns must contain the x- and y-coordinates with",
          " names 'xcoord' and 'ycoord'."
        )
      }
      if (any(colnames(covariate)[-(1:2)] == "")) {
        stop(
          "When spatial covariates are given as a matrix/data.frame, ",
          "all columns must be named."
        )
      }
      Z <- data.table::data.table(
        xcoord = X$x,
        ycoord = X$y,
        type_obs = as.integer(X$marks)
      )
      # q is the number of covariates including intercept. These maitrices are
      # not expected to include intercepts, but should contain the two spatial
      # coordinates of the points. This covariate contains q-1 coviariates and
      # 2 spatial coordinates (q+1 columns in total).
      q <- ncol(covariate) - 1
      covar_names <- colnames(covariate)[-(1:2)]
      spat_cov <- merge(spat_cov, covariate, by = c("xcoord", "ycoord"))
      rep_col <- rep(which(colnames(spat_cov) %in% covar_names), each = p - 1)

      rep_calc <- do.call(
        "cbind",
        replicate(q - 1, spat_cov[, 5:(p + 3)], simplify = F)
      )
      cov_mat <- spat_cov[, ..rep_col] * rep_calc
      colnames(cov_mat) <- paste(
        colnames(cov_mat),
        rep(1:(p - 1), q - 1),
        sep = ":"
      )
      spat_cov[, c(covar_names) := NULL]
      spat_cov <- cbind(spat_cov, cov_mat)
    }
    if (is.list(covariate) & !spatstat.geom::is.im(covariate) & !is.data.frame(covariate)) {
      check_im <- unlist(lapply(covariate, spatstat.geom::is.im))
      check_fct <- unlist(lapply(covariate, is.function))
      if (!all(check_im | check_fct)) {
        stop(
          "When coviarate is given as list, ",
          "all elements must be images or functions."
        )
      }
      q <- length(covariate) + 1
      for (i in 1:length(covariate)) {
        covar_name <- ifelse(
          is.null(names(covariate)[i]),
          paste("Covariate", i),
          names(covariate)[i]
        )
        if (spatstat.geom::is.im(covariate[[i]])) {
          spat_cov <- covar.im(X, covariate, covar_name, spat_cov)
        }
        if (is.function(covariate[[i]])) {
          spat_cov <- covar.function(X, covariate[[i]], covar_name, spat_cov)
        }
      }
    }
  }

  if (Poisson) {
    w <- spat_cov
  } else {
    Int_cov <- cbind(pts, Int_mat)
    w <- merge(spat_cov, Int_cov, by = c("xcoord", "ycoord", "type_obs", "j"))
  }

  w <- w[order(type_obs, xcoord, ycoord, j)]
  return(list(w = w, q = q))
}
