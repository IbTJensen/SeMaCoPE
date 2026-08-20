# Fast reimplementation of Covariate_setup() -----------------------------------
# Kept in its own file (rather than replacing the original in SemiMarkov.R) so
# the two versions can be benchmarked and validated against each other, e.g.
#   identical(Covariate_setup(X, Xis, nis, covariate, R_within, R_between, sat,
#                              Poisson, mark.pp)$w,
#             Covariate_setup_fast(X, Xis, nis, covariate, R_within, R_between,
#                                   sat, Poisson, mark.pp)$w)
#
# Five fixes on top of the original, none of which change the value, names,
# or order of any column of the returned `w` (only re-derived internally).
# Lettering (A)-(D) matches the original suggestions; the gc() removal wasn't
# part of that list and is lettered (E) here.
#
# (A) R/SemiMarkov.R:132-151 -- the (k, l, i) triple loop recomputes the full
#     "neighbours of type l within radius R[k,l]" sweep once per k, even
#     though R[k,l] only depends on k through the boolean (k == l): it is
#     R_within on the diagonal and R_between everywhere else. So for a given
#     l there are only 2 distinct radii to ever evaluate (within/between),
#     not p. Below, that sweep is computed once per (l, radius) and reused
#     for every k that shares the same radius, cutting the number of
#     closepairs()/crosspairs() calls from p^3 to 2*p^2. Columns are still
#     added to `Neighbours` in the exact original (k outer, l inner) order,
#     since downstream code (the Int_mat_list construction) indexes into
#     these columns positionally and assumes that order.
#
# (B) R/SemiMarkov.R:237-246 -- `DeltaS_maybe` was merged into `prelim_dt`
#     one (k, j) combination at a time (p^2 sequential merge.data.table()
#     calls), with `prelim_dt` getting wider -- and the join more expensive
#     -- on every iteration. Replaced with a single dcast() of
#     `Delta_S_maybe` into wide form (one row per (xcoord, ycoord, type_obs),
#     one column per "DeltaS_kj"), followed by exactly one merge into
#     `prelim_dt`. `setcolorder()` after the dcast pins the column order to
#     the exact k-outer/j-inner sequence the original loop produced, since
#     downstream code depends on that positional order (see (A)); dcast's
#     own column ordering is not relied upon.
#
# (C) R/SemiMarkov.R:253-264 -- `s <- as.matrix(prelim_dt[, ..m])` and
#     `Delta_S <- as.matrix(prelim_dt[, ..m])` do not depend on the loop
#     variable `j`, but were being re-extracted from `prelim_dt` on every one
#     of the p iterations. They're pulled out of the loop and computed once;
#     each iteration now only copies and zeroes the already-in-memory
#     matrices instead of re-subsetting/re-coercing the data.table each time.
#
# (D) R/SemiMarkov.R:70-73, inside DeltaS_init() -- `dummy_neighbours` was
#     built by calling apply(kl_comb, 1, ...) to construct one full n-row
#     copy of `Neighbours[, 1:3]` per row of `kl_comb` (up to p^2 of them),
#     purely to tag each copy with a different constant (Neighbour_type, l)
#     pair, then rbindlist()-ing all of those n-row tables back together.
#     Replaced with a single vectorized construction that repeats the row
#     indices of `Neighbours[, 1:3]` and the `kl_comb` columns directly,
#     producing the same rows in the same order in one allocation instead of
#     nrow(kl_comb) separate ones.
#
# (E) Profiling showed the explicit gc() calls scattered through DeltaS_init()
#     (R/SemiMarkov.R:40,51,53,55,57,88) and Covariate_setup() (line 168)
#     accounted for 54% of self-time on a p=8, n=1500 test case, dwarfing (A)
#     and (C) combined -- each gc() forces a full synchronous garbage
#     collection, and DeltaS_init() is called twice per Covariate_setup()
#     call, for 13 explicit collections per call. R's automatic GC already
#     runs when needed, so forcing it this often just adds fixed overhead
#     without reducing peak memory. DeltaS_init_fast() below is a verbatim
#     copy of DeltaS_init() with only the gc() calls removed, and the one in
#     Covariate_setup_fast() (after `rm(sat_mat)`) is removed the same way.
# Verbatim copy of DeltaS_init() (R/SemiMarkov.R:6-102) with the explicit
# gc() calls removed -- see note (E) above.
DeltaS_init_fast <- function(All_neighbours_between, Neighbours, p, between) {
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

  n_pts <- nrow(Neighbours)
  n_comb <- nrow(kl_comb)
  dummy_neighbours <- data.table::data.table(
    Neighbours[rep(seq_len(n_pts), times = n_comb), 1:3],
    Neighbour_type = rep(kl_comb$Neighbour_type, each = n_pts),
    l = rep(kl_comb$l, each = n_pts)
  )
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

Covariate_setup_fast <- function(
  X, Xis, nis, covariate, R_within,
  R_between, sat, Poisson, mark.pp
) {
  p <- length(Xis)
  n <- X$n
  sat_i <- rep(sat, p)
  sat_all <- rep(sat_i, nis)
  prelim_dt <- data.table::data.table(
    xcoord = X$x,
    ycoord = X$y,
    type_obs = as.integer(X$marks)
  )

  R <- matrix(R_between, p, p)
  diag(R) <- R_within

  if (!Poisson) {
    # Calculate the number of neighbours of each type each point has.
    Neighbours <- data.table::data.table(
      xcoord = X$x,
      ycoord = X$y,
      type_obs = as.integer(X$marks)
    )

    # For a fixed l, this is exactly the original inner (i in 1:p) sweep,
    # parameterised by the radius instead of implicitly via R[k, l].
    neighbours_of_type_l <- function(l, Xl, radius) {
      R_l_close_type_l_neighbours <- c()
      for (i in 1:p) {
        Xi <- Xis[[i]]
        if (i == l) {
          cc <- spatstat.geom::closepairs(Xi, rmax = radius)
        }
        if (i != l) {
          cc <- spatstat.geom::crosspairs(Xi, Xl, rmax = radius)
        }
        N <- table(factor(cc$i, levels = 1:nis[i]))
        R_l_close_type_l_neighbours <- c(R_l_close_type_l_neighbours, N)
      }
      R_l_close_type_l_neighbours
    }

    # Compute each of the (at most) 2 distinct radius outcomes per l once,
    # instead of once per (k, l) pair.
    N_within <- vector("list", p)
    N_between <- vector("list", p)
    for (l in 1:p) {
      Xl <- Xis[[l]]
      N_within[[l]] <- neighbours_of_type_l(l, Xl, R_within)
      N_between[[l]] <- neighbours_of_type_l(l, Xl, R_between)
    }

    for (k in 1:p) {
      for (l in 1:p) {
        R_kl_close_type_l_neighbours <- if (k == l) N_within[[l]] else N_between[[l]]
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

    DS_between <- DeltaS_init_fast(
      All_neighbours_between,
      Neighbours,
      p = p,
      between = T
    )
    DS_within <- DeltaS_init_fast(
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
    Delta_S_maybe[, kl_name := paste0("DeltaS_", Neighbour_type, l)]
    DeltaS_wide <- data.table::dcast(
      Delta_S_maybe,
      xcoord + ycoord + type_obs ~ kl_name,
      value.var = "DeltaS_kl"
    )
    kl_order <- paste0("DeltaS_", rep(1:p, each = p), rep(1:p, p))
    data.table::setcolorder(
      DeltaS_wide,
      c("xcoord", "ycoord", "type_obs", kl_order)
    )
    prelim_dt <- data.table::merge.data.table(
      prelim_dt,
      DeltaS_wide,
      by = c("xcoord", "ycoord", "type_obs")
    )
  }

  # Construct interaction matrix
  pts <- lapply(1:p, function(j) {
    data.table::data.table(prelim_dt[, 1:3], j = j)
  })
  pts <- data.table::rbindlist(pts)
  if (!Poisson) {
    m <- 3 + 1:p^2
    s_full <- as.matrix(prelim_dt[, ..m])
    m <- 3 + (p^2 + 1):(2 * p^2)
    Delta_S_full <- as.matrix(prelim_dt[, ..m])

    Int_mat_list <- list()
    for (j in 1:p) {
      s_idx <- 1:p + (j - 1) * p
      Delta_idx <- j + seq(1, p^2, p) - 1
      s <- s_full
      Delta_S <- Delta_S_full
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
