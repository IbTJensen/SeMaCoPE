#' @importFrom data.table data.table

# Calculate ΔS_ij for i,j=1,...,p ----------------------------------------------
DeltaS_init <- function(All_neighbours_between, Neighbours, p, between) {
  All_neighbours <- data.table::merge.data.table(
    x = Neighbours,
    y = All_neighbours_between,
    by = c("Neighbour_x", "Neighbour_y", "Neighbour_type")
  )
  data.table::setcolorder(All_neighbours, c("xcoord", "ycoord", "type_obs"))
  colnames(All_neighbours)[-(1:6)] <- paste(
    rep(1:p, each = p),
    rep(1:p, p),
    sep = "-"
  )

  within_idx <- 6 + seq(1, p^2, p) + 1:p - 1
  if (between) {
    between_idx <- (1:ncol(All_neighbours))[-within_idx]
    All_neighbours <- All_neighbours[, ..between_idx]
  }

  if (!between) {
    within_idx <- c(1:6, within_idx)
    All_neighbours <- All_neighbours[, ..within_idx]
  }

  strsplit_first_entry <- function(x, split, i) {
    unlist(lapply(strsplit(x, split = split), function(x) x[i]))
  }

  All_neighbours <- data.table::melt(
    data = All_neighbours,
    id.vars = 1:6,
    variable.name = "l",
    value.name = "s_Strauss_kl_v"
  )
  # Each row of All_neighbours contains points (u,i), its neighbour (v,k)
  # and s_kl(v, x_l\v) for k,l=1,...,p. Note that neighbour_type = k

  Int_types <- paste(rep(1:p, each = p), rep(1:p, p), sep = "-")
  Int_types1 <- strsplit_first_entry(Int_types, "-", 1)
  Int_types2 <- strsplit_first_entry(Int_types, "-", 2)
  names(Int_types1) <- Int_types
  names(Int_types2) <- Int_types

  All_neighbours[, l := as.character(l)]
  All_neighbours[, l2 := Int_types1[l]]
  All_neighbours <- All_neighbours[Neighbour_type == l2]
  All_neighbours[, ":="(l = as.integer(Int_types2[l]), l2 = NULL)]

  # For points with no neighbours, create a dummy point of each type,
  # where s_i is zero for all i
  kl_comb <- expand.grid(Neighbour_type = 1:p, l = 1:p)
  wtn_idx <- kl_comb$l == kl_comb$Neighbour_type
  if (between) {
    kl_comb <- kl_comb[!wtn_idx, ]
  }
  if (!between) {
    kl_comb <- kl_comb[wtn_idx, ]
  }

  apply(kl_comb, 1, function(x) {
    data.table::data.table(Neighbours[, 1:3], Neighbour_type = x[1], l = x[2])
  }) -> dummy_neighbours
  dummy_neighbours <- data.table::rbindlist(dummy_neighbours)
  colnames(dummy_neighbours) <- c(
    "xcoord",
    "ycoord",
    "type_obs",
    "Neighbour_type",
    "l"
  )

  All_neighbours <- data.table::merge.data.table(
    x = All_neighbours,
    y = dummy_neighbours,
    by = c("xcoord", "ycoord", "type_obs", "Neighbour_type", "l"),
    all = T
  )
  gc()
  All_neighbours[is.na(s_Strauss_kl_v), s_Strauss_kl_v := 0]
  data.table::setcolorder(
    All_neighbours,
    c(
      "xcoord",
      "ycoord",
      "type_obs",
      "Neighbour_x",
      "Neighbour_y",
      "Neighbour_type"
    )
  )
  return(All_neighbours)
}

# Sets up covariate data -------------------------------------------------------
Covariate_setup <- function(
  X, Xis, nis, covariate, R_within, 
  R_between, sat, Poisson, mark.pp
) {
  p <- length(Xis)
  n <- X$n
  # sat_i <- sat*nis/nis[p]
  sat_i <- rep(sat, p)
  sat_all <- rep(sat_i, nis)
  # The data.table prelim_dt will end up containing s_kl and Delta S_kl in each
  # point (i.e. the ingredients necessary to construct the w_i(u)'s).
  prelim_dt <- data.table::data.table(
    xcoord = X$x,
    ycoord = X$y,
    type_obs = as.integer(X$marks)
  )

  R <- matrix(R_between, p, p)
  diag(R) <- R_within

  if (!Poisson) {
    # Calculate the number of neighbours of each type each point has
    Neighbours <- data.table::data.table(
      xcoord = X$x,
      ycoord = X$y,
      type_obs = as.integer(X$marks)
    )
    for (k in 1:p) {
      for (l in 1:p) {
        Xl <- Xis[[l]]
        R_kl_close_type_l_neighbours <- c()
        for (i in 1:p) {
          Xi <- Xis[[i]]
          if (i == l) {
            cc <- spatstat.geom::closepairs(Xi, rmax = R[k, l])
          }
          if (i != l) {
            cc <- spatstat.geom::crosspairs(Xi, Xl, rmax = R[k, l])
          }
          N <- table(factor(cc$i, levels = 1:nis[i]))
          R_kl_close_type_l_neighbours <- c(R_kl_close_type_l_neighbours, N)
        }
        Neighbours[, a := R_kl_close_type_l_neighbours]
        m <- ncol(Neighbours)
        colnames(Neighbours)[m] <- paste0("N_", k, l)
      }
    }
    sat_mat <- matrix(rep(sat_i, p), nrow = n, ncol = p^2, byrow = T)
    colnames(sat_mat) <- paste(
      "sat",
      paste0(rep(1:p, each = p), rep(1:p, p)),
      sep = "_"
    )

    # Calculate s_ik for all combinations. Note that we assume that
    # R_kl = R_lk, and that s_kl(u) denotes the number of R_kl-close type l
    # points of u (up to saturation). Thus s_kl = s_il for all k,i. Below,
    # s_i is such that s_i = s_ki for all k.
    s_mat <- data.table::data.table(Neighbours[, -(1:3)])
    colnames(s_mat) <- gsub("N", "s", colnames(s_mat))
    s_mat[s_mat > sat_mat] <- sat_mat[s_mat > sat_mat]
    prelim_dt <- cbind(prelim_dt, s_mat)
    rm(sat_mat)
    gc()

    # Constructing all pairs of within and between neighbours
    All_neighbours_between_kl <- list()
    All_neighbours_within_kk <- list()
    j <- 1
    for (k in 1:p) {
      Xk <- Xis[[k]]
      for (l in 1:p) {
        Xl <- Xis[[l]]
        if (l == k) {
          cc_between <- spatstat.geom::closepairs(Xl, rmax = R_between)
          cc_within <- spatstat.geom::closepairs(Xl, rmax = R_within)
        }
        if (l != k) {
          cc_between <- spatstat.geom::crosspairs(Xk, Xl, rmax = R_between)
          cc_within <- spatstat.geom::crosspairs(Xk, Xl, rmax = R_within)
        }
        dt <- data.table::data.table(
          xcoord = cc_between$xi,
          ycoord = cc_between$yi,
          type_obs = k,
          Neighbour_x = cc_between$xj,
          Neighbour_y = cc_between$yj,
          Neighbour_type = l
        )
        All_neighbours_between_kl[[j]] <- dt

        dt <- data.table::data.table(
          xcoord = cc_within$xi,
          ycoord = cc_within$yi,
          type_obs = k,
          Neighbour_x = cc_within$xj,
          Neighbour_y = cc_within$yj,
          Neighbour_type = l
        )
        All_neighbours_within_kk[[j]] <- dt
        j <- j + 1
      }
    }

    All_neighbours_between <- data.table::rbindlist(All_neighbours_between_kl)
    All_neighbours_within <- data.table::rbindlist(All_neighbours_within_kk)

    colnames(Neighbours)[1:3] <- c(
      "Neighbour_x",
      "Neighbour_y",
      "Neighbour_type"
    )

    DS_between <- DeltaS_init(
      All_neighbours_between,
      Neighbours,
      p = p,
      between = T
    )
    DS_within <- DeltaS_init(
      All_neighbours_within,
      Neighbours,
      p = p,
      between = F
    )
    DS <- rbind(DS_between, DS_within)
    DS[, l := as.numeric(as.character(l))]
    DS[, sat := sat_i[as.numeric(l)]]
    DS[,
      .(DeltaS_kl = sum(!is.na(Neighbour_x) & s_Strauss_kl_v <= sat - 1)),
      by = list(xcoord, ycoord, type_obs, Neighbour_type, l)
    ] -> Delta_S_maybe
    for (k in 1:p) {
      for (j in 1:p) {
        data.table::merge.data.table(
          prelim_dt,
          Delta_S_maybe[Neighbour_type == k & l == j, -(4:5)],
          by = c("xcoord", "ycoord", "type_obs")
        ) -> prelim_dt
        colnames(prelim_dt)[ncol(prelim_dt)] <- paste0("DeltaS_", k, j)
      }
    }
  }

  # Construct interaction matrix
  pts <- lapply(1:p, function(j) data.table::data.table(prelim_dt[, 1:3], j = j))
  pts <- data.table::rbindlist(pts)
  if (!Poisson) {
    Int_mat_list <- list()
    for (j in 1:p) {
      m <- 3 + 1:p^2
      s <- as.matrix(prelim_dt[, ..m])
      m <- 3 + (p^2 + 1):(2 * p^2)
      Delta_S <- as.matrix(prelim_dt[, ..m])
      s_idx <- 1:p + (j - 1) * p
      Delta_idx <- j + seq(1, p^2, p) - 1
      s[, -s_idx] <- 0
      Delta_S[, -Delta_idx] <- 0
      Int_mat_list[[j]] <- s + Delta_S
    }
    Int_mat <- Reduce(rbind, Int_mat_list)
    colnames(Int_mat) <- paste(rep(1:p, each = p), rep(1:p, p), sep = "-")
  }

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
      # rep_col <- c(rep(3:(q+1), each = p-1))
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

beta_estimation <- function(w, betastart, p, q, Poisson) {
  Wmat <- as.matrix(w[, -(1:4), with = FALSE])

  # Group id per row (one group per observed point, i.e. per (xcoord, ycoord,
  # type_obs) combination), computed once instead of re-derived every call.
  grp_dt <- w[, .(xcoord, ycoord, type_obs)]
  grp_dt[, grp := .GRP, by = .(xcoord, ycoord, type_obs)]
  grp <- grp_dt$grp

  # Row indices of the observed ("true") type, and their group ids, also
  # computed once.
  true_idx <- which(w$type_obs == w$j)
  true_grp <- grp[true_idx]

  loglik <- function(beta) {
    cond_intensity <- as.vector(exp(Wmat %*% beta))
    Lambda <- rowsum(cond_intensity, grp)[, 1]
    prob_true <- cond_intensity[true_idx] / Lambda[true_grp]
    sum(log(prob_true[!is.na(prob_true)]))
  }

  score <- function(beta){
    cond_intensity <- as.vector(exp(Wmat %*% beta))
    Lambda <- rowsum(cond_intensity, grp)[, 1]
    prob <- cond_intensity / Lambda[grp]
    hmat <- Wmat[true_idx,] - rowsum(Wmat*prob, grp)
    colSums(hmat, na.rm = TRUE)
  }

  obj <- optim(
    par = betastart,
    fn = loglik,
    gr = score,
    method = "BFGS",
    control = list(fnscale = -1)
  )

  return(obj)
}

# beta_estimation <- function(w, betastart, p, q, Poisson) {
#   loglik <- function(beta) {
#     pred <- type_pred(w, beta)
#     pred <- pred[type_obs == j]
#     loglik <- sum(log(pred$prob[!is.na(pred$prob)]))
#     return(loglik)
#   }

#   obj <- optim(
#     par = betastart,
#     fn = loglik,
#     method = "BFGS",
#     control = list(fnscale = -1)
#   )
#   return(obj)
# }

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
  R_between, sat = Inf, standardize = TRUE, Poisson = FALSE
) {
  if (!(is.integer(X$marks) | is.factor(X$marks))) {
    stop("Marks of point process must be of type either integer or factor.")
  }

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
  p <- length(mark.pp)
  Xis = list()
  nis = rep(0, p)
  for (i in 1:p) {
    Xis[[i]] = X[mark.pp[i] == X$marks]
    nis[i] = Xis[[i]]$n
  }

  # w <- Covariate_setup_fast(
  #   X, Xis, nis, covariate, R_within, 
  #   R_between, sat, Poisson, mark.pp
  # )
  w <- Covariate_setup_rcpp(
    X, Xis, nis, covariate, R_within, 
    R_between, sat, Poisson, mark.pp
  )
  q <- w$q
  w <- w$w

  # a <- Sys.time()
  # if (!is.null(edgecorrection)) {
  #   erodedwindow = spatstat.geom::erosion(X$window, edgecorrection)
  #   pts_in_window <- spatstat.geom::inside.owin(x = w$xcoord, y = w$ycoord, w = erodedwindow)
  #   w <- w[pts_in_window]
  #   X <- X[spatstat.geom::inside.owin(X, w = erodedwindow), ]
  #   for (i in 1:p) {
  #     pts_in_window <- spatstat.geom::inside.owin(Xis[[i]], w = erodedwindow)
  #     Xis[[i]] = Xis[[i]][pts_in_window, ]
  #     nis[i] = sum(pts_in_window)
  #   }
  # }
  # Sys.time() - a

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
    offsets <- rep(0, q * (p - 1) + p * (p + 1) / 2)
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
  opt <- beta_estimation_rcpp(w, betastart, p, q, Poisson)
  betahat <- opt$par
  # w <- w_raw
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
    sat = sat,
    Xis = Xis,
    w = w_raw,
    q = q,
    maximum_log_likelihood = opt$value,
    pred = pred,
    pred_no_int = pred_no_int
  )

  return(out)
}

# Estimate standard errors -----------------------------------------------------
Standard_error_matrix <- function(X, w, betahat, R_within, R_between, sat) {
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
  Int_range <- ifelse(sat == Inf, Int_range, 2 * Int_range)
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

  # cc <- spatstat.geom::closepairs(X_new, rmax = Int_range)
  # gc()
  # Hu <- H[cc$i, ]
  # Hv <- H[cc$j, ]
  # Sigma_term <- t(Hu) %*% Hv
  # rm(Hu, Hv)
  # gc()

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
#' @param sat Indicates the saturation parameter. If set to Inf, a Strauss
#' model will be fitted. If a single value is provided, this value will be used.
#' If a vector is used, a grid search is carried out, and the value that
#' maximises the log composite likelihood is used.
#' @param Poisson If TRUE, a Poisson process is fitted. If FALSE interaction
#' terms will be included.
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
  X, covariate, edgecorrection = NULL, R_within,
  R_between, sat = Inf, standardize = TRUE, Poisson = FALSE
) {
  if (length(R_within) == 1 & length(R_between) == 1 & length(sat) == 1) {
    a <- Sys.time()
    S <- SemiMarkov_fixed_R(
      X,
      covariate = covariate,
      edgecorrection = edgecorrection,
      R_within = R_within,
      R_between = R_between,
      sat = sat,
      standardize = standardize,
      Poisson = Poisson
    )
    Sys.time() - a

    a <- Sys.time()
    std_err <- Standard_error_matrix(
      X = X,
      w = S$w,
      betahat = S$betahat,
      R_within = R_within,
      R_between = R_between,
      sat = sat
    )
    Sys.time() - a

    out <- c(S, std_err)
    return(out)
  }

  opt <- expand.grid(R_within = R_within, R_between = R_between, sat = sat)
  R_omit <- opt$R_within == opt$R_between
  opt <- opt[!R_omit, ]
  res_list <- list()
  likelihoods <- rep(NA, nrow(opt))
  for (i in 1:nrow(opt)) {
    cat(paste("Fitting model", i, "out of", nrow(opt)), "\r")
    S <- SemiMarkov_fixed_R(
      X,
      covariate = covariate,
      edgecorrection = edgecorrection,
      R_within = opt[i, "R_within"],
      R_between = opt[i, "R_between"],
      sat = opt[i, "sat"],
      standardize = standardize,
      Poisson = Poisson
    )
    res_list[[i]] <- S
    likelihoods[i] <- S$maximum_log_likelihood
  }

  w <- which.max(likelihoods[likelihoods != 0])
  S <- res_list[[w]]
  std_err <- Standard_error_matrix(
    X,
    S$w,
    S$betahat,
    opt[w, "R_within"],
    opt[w, "R_between"],
    opt[w, "sat"]
  )
  out <- c(S, std_err)
  return(out)
}

# Check second order term ------------------------------------------------------
# compute the sensitivity and second order term in the true value. This is
# for model evaluation purposes
Covariance_true_val <- function(
  X, covariate, edgecorrection = NULL, R_within,
  R_between, sat = Inf, Poisson = FALSE, true.param
) {
  if (!(is.integer(X$marks) | is.factor(X$marks))) {
    stop("Marks of point process must be of type either integer or factor.")
  }

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
  p <- length(mark.pp)
  Xis = list()
  nis = rep(0, p)
  for (i in 1:p) {
    Xis[[i]] = X[mark.pp[i] == X$marks]
    nis[i] = Xis[[i]]$n
  }

  w <- Covariate_setup(
    X,
    Xis,
    nis,
    covariate,
    R_within,
    R_between,
    sat,
    Poisson,
    mark.pp
  )
  q <- w$q
  w <- w$w

  if (!is.null(edgecorrection)) {
    erodedwindow = spatstat.geom::erosion(X$window, edgecorrection)
    pts_in_window <- spatstat.geom::inside.owin(x = w$xcoord, y = w$ycoord, w = erodedwindow)
    w <- w[pts_in_window]
    X <- X[spatstat.geom::inside.owin(X, w = erodedwindow), ]
    for (i in 1:p) {
      pts_in_window <- spatstat.geom::inside.owin(Xis[[i]], w = erodedwindow)
      Xis[[i]] = Xis[[i]][pts_in_window, ]
      nis[i] = sum(pts_in_window)
    }
  }

  w_raw <- data.table::data.table(w)
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

  var_cols <- 5:ncol(w)
  var_col_names <- colnames(w)[var_cols]

  # Fitted log-lambda (up to non-parametric factor)
  log_lambda <- as.matrix(w[, -(1:4)]) %*% true.param
  w[, lambda := exp(log_lambda)]
  w_Lambda <- w[, .(Lambda = sum(lambda)), list(xcoord, ycoord, type_obs)]
  w <- merge(w, w_Lambda, by = c("xcoord", "ycoord", "type_obs"))
  w[, prob := lambda / Lambda]

  h <- w[,
    lapply(.SD, function(x) x - sum(x * prob)),
    by = list(xcoord, ycoord, type_obs),
    .SDcols = var_col_names
  ]
  h <- cbind(w[, 1:4], h[, -(1:3)])

  # Matrix-matrix product
  w_mat <- as.matrix(w[, ..var_cols])
  h_mat <- as.matrix(h[, ..var_cols])
  probs <- w$prob
  S <- t(w_mat * probs) %*% h_mat

  # Second order term
  Int_range <- max(R_within, R_between)
  Int_range <- ifelse(sat == Inf, Int_range, 2 * Int_range)
  h_type_obs <- h[type_obs == j]
  H <- as.matrix(h_type_obs[match(X$x, h_type_obs$xcoord), -(1:4)])
  cc <- spatstat.geom::closepairs(X, rmax = Int_range)
  Hu <- H[cc$i, ]
  Hv <- H[cc$j, ]
  Sigma_term <- t(Hu) %*% Hv

  rownames(Sigma_term) <- names(betahat)
  colnames(Sigma_term) <- names(betahat)
  rownames(S) <- names(betahat)
  colnames(S) <- names(betahat)

  out <- list(S = S, Sigma_term = Sigma_term)
  return(out)
}
