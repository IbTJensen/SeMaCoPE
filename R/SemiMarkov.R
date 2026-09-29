#' @importFrom data.table data.table

type_pred <- function(w, beta) {
  cond_intensity <- exp(as.matrix(w[, -(1:4), with = F]) %*% beta)
  pts_lam <- w[, 1:4]
  pts_lam[, lambda := cond_intensity]
  Lam <- pts_lam[, .(Lambda = sum(lambda)), list(xcoord, ycoord, type_obs)]
  pts_lam <- data.table::merge.data.table(
    pts_lam,
    Lam,
    by = c("xcoord", "ycoord", "type_obs")
  )
  pts_lam[, prob := lambda / Lambda]
  pts_lam <- pts_lam[order(type_obs, xcoord, ycoord, j)]
  return(pts_lam)
}

# Fitting Semi parametric Markov model for fixed R and c -----------------------
SemiMarkov_fixed_R <- function(
  X, covariate, edgecorrection = NULL, R_within,
  R_between, sat_within = Inf, sat_between = Inf,
  standardize = TRUE, Poisson = FALSE, quiet = FALSE
) {
  if (!is.integer(X$marks) & !is.factor(X$marks) & !is.character(X$marks)) {
    stop("Marks of X must characters, integers or factor.")
  }

  p <- length(unique(X$marks))
  if (is.integer(X$marks)) {
    marks <- sort(unique(X$marks))
    X$marks <- factor(X$marks, levels = marks[c(p, 1:(p - 1))])
  }

  if(is.character(X$marks)){
    X$marks <- factor(X$marks, levels = sort(unique(X$marks)))
  }

  X$marks <- factor(X$marks, levels = levels(X$marks)[c(2:p, 1)])

  # Jitters duplicated points by a small distance
  m <- min(R_within, R_between)
  duplicates <- duplicated(data.frame(X$x, X$y))
  d <- sum(duplicates)

  perturb_x <- runif(d, -m / 1000, m / 1000)
  perturb_y <- runif(d, -m / 1000, m / 1000)

  if (is.matrix(covariate) | is.data.frame(covariate)) {
    pts_X <- paste(X$x, X$y, sep = "-")
    pts_covariates <- paste(covariate$xcoord, covariate$ycoord, sep = "-")
    reorder <- match(pts_X, pts_covariates)
    covariate <- covariate[reorder, ]
    covariate$xcoord[duplicates] <- covariate$xcoord[duplicates] + perturb_x
    covariate$ycoord[duplicates] <- covariate$ycoord[duplicates] + perturb_y
  }

  X$x[duplicates] <- X$x[duplicates] + perturb_x
  X$y[duplicates] <- X$y[duplicates] + perturb_y

  mark.pp <- sort(unique(X$marks))
  Xis = list()
  nis = rep(0, p)
  for (i in 1:p) {
    Xis[[i]] = X[mark.pp[i] == X$marks]
    nis[i] = Xis[[i]]$n
  }

  w <- Covariate_setup_sat2_rcpp(
    X, Xis, nis, covariate, R_within, 
    R_between, sat_within, sat_between, Poisson, mark.pp
  )
  q <- w$q
  w <- w$w

  if (!is.null(edgecorrection)) {
    ec <- apply_edge_correction(X, Xis, nis, w, p, mark.pp, edgecorrection)
    X <- ec$X
    Xis <- ec$Xis
    nis <- ec$nis
    w <- ec$w
  }

  if (!Poisson) {
    Int_var_dt <- w[, (4 + q * (p - 1) + 1):ncol(w)]
    nms <- colnames(Int_var_dt)
    lapply(strsplit(nms, "-"), function(x) {
      if (x[1] > x[2]) {
        return(paste(x[2:1], collapse = "-"))
      } else {
        return(paste(x, collapse = "-"))
      }
    }) -> nms
    nms <- unlist(nms)
    w_sym <- data.table::data.table(t(rowsum(t(Int_var_dt), nms)))
    w <- cbind(w[, 1:(4 + q * (p - 1))], w_sym)
  }
  w_raw <- data.table::data.table(w)

  w_var <- w[, (4 + p):ncol(w)]
  if (standardize) {
    offsets <- colMeans(w_var, na.rm = T)
    scalings <- apply(w_var, 2, sd, na.rm = T)
  } else {
    offsets <- rep(0, (q - 1) * (p - 1) + p * (p + 1) / 2)
    scalings <- rep(1, (q - 1) * (p - 1) + p * (p + 1) / 2)
  }

  offset_mat <- matrix(
    offsets,
    nrow = nrow(w),
    ncol = length(offsets),
    byrow = T
  )
  scaling_mat <- matrix(
    scalings,
    nrow = nrow(w),
    ncol = length(offsets),
    byrow = T
  )

  w[, (4 + p):ncol(w)] <- (w_var - offset_mat) / scaling_mat
  back_transform <- c(rep(1, p - 1), scalings)

  response <- rep(1:p, nis)
  response <- factor(response)

  if (is.null(covariate)) {
    fit = VGAM::vglm(response ~ 1, family = VGAM::multinomial)
    betafitz = coef(fit)
  } else {
    covar_cols <- colnames(w)[
      grepl(":1", colnames(w)) &
        colnames(w) != "Intercept:1"
    ]
    Z <- as.matrix(w[order(j)][j == 1, ..covar_cols])
    colnames(Z) <- gsub(":1", "", colnames(Z))

    fit = VGAM::vglm(response ~ Z, family = VGAM::multinomial)
    betafitz = coef(fit)

    if (q == 2) {
      names(betafitz) <- gsub("Z", "Covariate", names(betafitz))
    }
    if (q > 2) {
      names(betafitz) <- gsub("Z", "", names(betafitz))
    }
  }
  predictions = predict(fit, type = "response") #predicted probabilities

  # Add zeros to interaction parameters as starting values
  if (!Poisson) {
    betastart <- c(betafitz, rep(0, (p + 1) * p / 2))
    names(betastart)[-(1:(q * (p - 1)))] <- paste(
      rep(1:p, p:1),
      sequence(p:1, from = 1:p),
      sep = "-"
    )
  } else {
    betastart <- betafitz
  }

  if (p == 2) {
    names(betastart)[1:q] <- paste(names(betastart)[1:q], "1", sep = ":")
  }

  # Estimate parameters
  opt <- beta_estimation_rcpp(w, betastart, p, q, Poisson, quiet = quiet)
  betahat <- opt$par
  w <- w_raw
  var_cols <- 5:ncol(w)
  var_col_names <- colnames(w)[var_cols]
  betahat <- betahat / back_transform
  pred <- type_pred(w, betahat)
  colnames(pred)[4] <- "type_pred"
  pred[, ":="(lambda = NULL, Lambda = NULL)]

  betahat_no_int <- c(betahat)
  betahat_no_int[-(1:(q * (p - 1)))] <- 0
  cond_no_int <- exp(as.matrix(w[, -c(1:4), with = F]) %*% betahat_no_int)

  pred_no_int <- w[, 1:4]
  pred_no_int[, lambda := cond_no_int]
  Lam <- pred_no_int[, .(Lambda = sum(lambda)), list(xcoord, ycoord, type_obs)]
  pred_no_int <- data.table::merge.data.table(
    pred_no_int,
    Lam,
    by = c("xcoord", "ycoord", "type_obs")
  )
  pred_no_int[, prob := lambda / Lambda]
  pred_no_int <- pred_no_int[order(type_obs, xcoord, ycoord, j)]

  colnames(pred_no_int)[4] <- "type_pred"
  pred_no_int[, ":="(lambda = NULL, Lambda = NULL)]

  out <- list(
    betahat = betahat,
    betahatpois = betafitz,
    predictions = predictions,
    ref_type = levels(X$marks)[1],
    R_within = R_within,
    R_between = R_between,
    sat_within = sat_within,
    sat_between = sat_between,
    Xis = Xis,
    w = w_raw,
    q = q,
    maximum_log_likelihood = opt$value,
    convergence = opt$convergence,
    convergence_message = opt$message,
    pred = pred,
    pred_no_int = pred_no_int
  )

  return(out)
}

# Estimate standard errors -----------------------------------------------------
Standard_error_matrix <- function(X, w, betahat, R_within, R_between, sat_within, sat_between) {
  # Fitted log-lambda (up to non-parametric factor)
  log_lambda <- as.matrix(w[, -(1:4)]) %*% betahat
  w[, lambda := exp(log_lambda)]
  w_Lambda <- w[, .(Lambda = sum(lambda)), list(xcoord, ycoord, type_obs)]
  w <- merge(w, w_Lambda, by = c("xcoord", "ycoord", "type_obs"))
  w[, prob := lambda / Lambda]

  var_cols <- 5:(ncol(w) - 3)
  var_col_names <- colnames(w)[var_cols]

  h <- w[,
    lapply(.SD, function(x) x - sum(x * prob)),
    by = list(xcoord, ycoord, type_obs),
    .SDcols = var_col_names
  ]
  h <- cbind(w[, 1:4], h[, -(1:3)])

  # Sensitivity matrix
  w_mat <- as.matrix(na.omit(w)[, ..var_cols])
  h_mat <- as.matrix(na.omit(h)[, ..var_cols])
  probs <- na.omit(w)$prob
  S <- t(w_mat * probs) %*% h_mat

  # Second order term
  Int_range <- max(R_within, R_between)
  Int_range <- ifelse(sat_within == Inf & sat_between == Inf, Int_range, 2 * Int_range)
  h_type_obs <- na.omit(h)[type_obs == j]

  X_new <- spatstat.geom::ppp(
    x = h_type_obs$xcoord,
    y = h_type_obs$ycoord,
    window = X$window,
    marks = factor(h_type_obs$type_obs)
  )

  H <- as.matrix(h_type_obs[, -(1:4)])
  n <- nrow(H)
  cc <- spatstat.geom::closepairs(X_new, rmax = Int_range)
  W <- Matrix::sparseMatrix(i = cc$i, j = cc$j, x = 1, dims = c(n, n))
  Sigma_term <- as.matrix(t(H) %*% (W %*% H))

  rownames(Sigma_term) <- names(betahat)
  colnames(Sigma_term) <- names(betahat)
  rownames(S) <- names(betahat)
  colnames(S) <- names(betahat)

  Sigma <- S + Sigma_term
  S_inv <- solve(S)
  V <- S_inv %*% Sigma %*% S_inv
  std_err <- sqrt(diag(V))

  CI <- data.table::data.table(
    Covariate = names(betahat),
    Estimate = betahat,
    Lower_CI = betahat - 1.96 * std_err,
    Upper_CI = betahat + 1.96 * std_err
  )

  z_stats <- betahat / std_err
  p_vals <- 2 * pnorm(-abs(z_stats))

  out <- list(
    std_err = std_err,
    std_err_S = sqrt(diag(S_inv)),
    Second_order_term = Sigma_term,
    Sensitivity = S,
    CI = CI,
    p_vals = p_vals,
    V = V,
    h = h_type_obs
  )

  return(out)
}

# Semi-parametric Markov model -------------------------------------------------
loglik_internal <- function(
  X, covariate, edgecorrection, R_within, 
  R_between, sat_within, sat_between, standardize, Poisson
) {
  S <- SemiMarkov_fixed_R(
    X,
    covariate,
    edgecorrection,
    R_within,
    R_between,
    sat_within,
    sat_between,
    standardize,
    Poisson,
    quiet = TRUE
  )

  return(S$maximum_log_likelihood)
}


#' Fits semi-parametric Markov model to multitype point pattern data with a
#' fixed interaction radius and saturation parameter.
#' @param X Multiyupe point pattern data. Must be a spatstat ppp object with
#' the point type as marks.
#' @param covariate Covariate(s) to be included in the model. Must be either
#' a matrix/data.frame (in the form of a design matrix),a spatstat pixel image
#' object, or a function that takes an x- and y-coordinates as inputs. If NULL
#' is provided, only an intercept will be fitted.
#' @param edgecorrection Indicates the distance with which the observation
#' windows should be eroded to account for edge effects. If NULL is provided,
#' no edge correction is performed.
#' @param R_within The interaction range within types. If a single value is
#' provided, this value will be used. If a vector is used, a grid search is
#' carried out, and the value that maximises the log composite likelihood is
#' used.
#' @param R_between The interaction range between types. If a single value is
#' provided, this value will be used. If a vector is used, a grid search is
#' carried out, and the value that maximises the log composite likelihood is
#' used.
#' @param standardize Logical argumenet indicating whether or not
#' standardization should be performed internally on the covariates. For
#' interpretability, parameter estimates are transformed back to take this
#' into account. It is recommended to leave as TRUE.
#' @param sat Indicates the common saturation parameter among all types.
#' For seperate within-type and between-type saturation parameters, set sat = NULL.
#' If set to Inf, a Strauss model will be fitted. If a single value is provided, 
#' this value will be used. If a vector is used, a grid search is carried out, 
#' and the value that maximises the log composite likelihood is used.
#' @param sat_within Indicates the within-type saturation parameter. 
#' To use this (and sat_between), set sat = NULL. If set to Inf, a Strauss model 
#' will be fitted. If a single value is provided, this value will be used.
#' If a vector is used, a grid search is carried out, and the value that
#' maximises the log composite likelihood is used.
#' @param sat_within Indicates the between-type saturation parameter. 
#' To use this (and sat_within), set sat = NULL. If set to Inf, a Strauss model 
#' will be fitted. If a single value is provided, this value will be used.
#' If a vector is used, a grid search is carried out, and the value that
#' maximises the log composite likelihood is used.
#' @param Poisson If TRUE, a Poisson process is fitted. If FALSE interaction
#' terms will be included.
#' @param quiet Should the function keep quiet about points with missing values
#' in the covariate being dropped?
#' @return A list that includes betahat (the parameter estimate of beta),
#' converg (convergence information passed from optim), ref_type (the reference
#' type (currently not implemented)), R_within (the interaction range between
#' types used by the final model), R_between (the interaction range between
#' types used by the final model), sat (the saturation parameter used by the
#' final model), Xis (a list of point processes of the different types), w (a
#' data table containing the w_i(u) vectors for each point, u in the point
#' processes), q (the number of covariates used in the model aside from
#' intercept), and maximum_log_likelihood (the maximum log composite
#' likelihood).
#' @export
SemiMarkov <- function(
  X, covariate, edgecorrection = NULL, R_within, R_between, 
  sat = Inf, sat_within = NULL, sat_between = NULL, standardize = TRUE, 
  Poisson = FALSE, ncores = 1, quiet = FALSE
) {
  if (!is.null(sat) & (!is.null(sat_within) | !is.null(sat_between))) {
    err_msg <- paste(
      "sat_within and sat_between must be NULL when sat is non-NULL."
    )
    stop(err_msg)
  }

  if (is.null(sat) & (is.null(sat_within) | is.null(sat_between))) {
    err_msg <- paste(
      "When sat = NULL both sat_within and sat_between must be non-NUll."
    )
    stop(err_msg)
  }

  if (!is.null(sat)) {
    sat_within <- sat
    sat_between <- sat
  }

  is_R_len1 <- length(R_within) == 1 & length(R_between) == 1
  is_sats_len1 <- length(sat_within) == 1 & length(sat_between == 1)

  if (is_R_len1 & is_sats_len1) {

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

  if (!is.null(sat)) {
    sat_omit <- opt$sat_within != opt$sat_between
    opt <- opt[!sat_omit, ]
  }

  if (ncores == 1) {
    likelihoods <- rep(NA, nrow(opt))
    for (i in 1:nrow(opt)) {
      cat(paste("Fitting model", i, "out of", nrow(opt)), "\r")
      likelihoods[i] <- loglik_internal(
        X,
        covariate = covariate,
        edgecorrection = edgecorrection,
        R_within = opt[i, "R_within"],
        R_between = opt[i, "R_between"],
        sat_within = opt[i, "sat_within"],
        sat_between = opt[i, "sat_between"],
        standardize = standardize,
        Poisson = Poisson
      )
    }
  }

  if (ncores > 1) {
    helper_fct <- function(i) {
      l <- loglik_internal(
        X,
        covariate = covariate,
        edgecorrection = edgecorrection,
        R_within = opt[i, "R_within"],
        R_between = opt[i, "R_between"],
        sat_within = opt[i, "sat_within"],
        sat_between = opt[i, "sat_between"],
        standardize = standardize,
        Poisson = Poisson
      )
      return(l)
    }

    future::plan(future::multisession, workers = ncores)
    lik_list <- furrr::future_map(
      1:nrow(opt),
      helper_fct,
      .progress = TRUE
    )

    likelihoods <- unlist(lik_list)
  }
  liks <- data.table::data.table(opt, likelihoods = likelihoods)

  w <- which.max(likelihoods[likelihoods != 0])
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

