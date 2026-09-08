# Compiled (Rcpp) reimplementation of Covariate_setup() ------------------------
# Kept in its own file, alongside Covariate_setup_fast() (the pure-R version
# with fixes A-E), so the versions can be benchmarked and validated against
# each other, e.g.
#   all.equal(Covariate_setup(X, Xis, nis, covariate, R_within, R_between,
#                             sat, Poisson, mark.pp)$w,
#             Covariate_setup_rcpp(X, Xis, nis, covariate, R_within, R_between,
#                                  sat, Poisson, mark.pp)$w)
#
# Same signature and same return value -- list(w = <data.table>, q = <int>) --
# as Covariate_setup(), so it is a drop-in replacement inside
# SemiMarkov_fixed_R().
#
# WHAT IS COMPILED AND WHAT IS NOT
# --------------------------------
# The spatial/combinatorial core -- neighbour counts, saturation, the DeltaS
# computation (which in the R version is DeltaS_init()'s merge/melt/join
# pipeline over every neighbour pair), and the Int_mat construction -- runs in
# src/covariate_setup_rcpp.cpp, compiled into the package DLL at install time
# and called through the generated RcppExports wrapper. See that file's header
# for the identity that collapses DeltaS_init() into two weighted radius-query
# passes.
#
# The covariate handling below (spatstat `im` lookup, covariate functions,
# user-supplied data.frames, the intercept dummies and their interactions)
# stays in R, unchanged from Covariate_setup(). It is cheap relative to the
# spatial work, and porting it would mean calling back into R regardless.
#
# IMPORTANT -- A PRE-EXISTING BUG THIS DOES NOT REPRODUCE
# -------------------------------------------------------
# Covariate_setup() (R/SemiMarkov.R:132-151) builds each neighbour-count
# column by concatenating over types:
#
#     for (i in 1:p) { ...; R_kl_close_type_l_neighbours <- c(R_kl..., N) }
#     Neighbours[, a := R_kl_close_type_l_neighbours]
#
# The vector is therefore in type-sorted order, while `Neighbours` rows are in
# the original order of X. Unless X's marks are already sorted, every count is
# attached to the wrong point, and the error propagates into s, DeltaS and
# Int_mat. Verified against a brute-force count: the R version matches only
# when X is sorted by mark.
#
# This implementation indexes points directly and is correct for any input
# ordering. Consequences:
#   * X sorted by mark    -> output is identical to Covariate_setup().
#   * X not sorted by mark -> output differs, and this one is the correct one.
# To reproduce the old (incorrect) numbers for comparison, sort X by mark
# before calling either version; both then agree.
Covariate_setup_rcpp <- function(
  X, Xis, nis, covariate, R_within,
  R_between, sat, Poisson, mark.pp
) {
  p <- length(Xis)
  n <- X$n
  if (length(sat) %notin% c(1, p)) {
    err_msg <- paste(
      "sat should be either a single number or a",
      "vector with one entry per point type."
    )
    stop(err_msg)
  }
  
  sat_i <- sat
  if(length(sat) == 1){
    sat_i <- rep(sat, p)
  }

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
      stop("Covariate_setup_rcpp(): every mark of X must appear in mark.pp.")
    }

    core <- covariate_setup_core_cpp(
      x = as.numeric(X$x),
      y = as.numeric(X$y),
      type = as.integer(type_code),
      p = p,
      R_within = R_within,
      R_between = R_between,
      sat_i = as.numeric(sat_i),
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
