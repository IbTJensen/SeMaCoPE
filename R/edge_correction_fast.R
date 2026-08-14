# Faster edge correction -------------------------------------------------------
# Drop-in replacement for the edge-correction block in SemiMarkov_fixed_R()
# (R/SemiMarkov.R:542-552). Kept in its own file so the two can be compared:
#
#   old <- { erodedwindow = erosion(X$window, edgecorrection); ... }
#   new <- apply_edge_correction(X, Xis, nis, w, p, mark.pp, edgecorrection)
#
# The original block performs n * p + 2n point-in-window tests:
#
#   inside.owin(x = w$xcoord, y = w$ycoord, w = erodedwindow)   # n * p rows
#   inside.owin(X, w = erodedwindow)                            # n points
#   for (i in 1:p) inside.owin(Xis[[i]], w = erodedwindow)      # n points total
#
# Two of those three are redundant:
#
# (1) `w` holds each observed point exactly p times -- once per candidate type
#     j -- so testing every row tests each distinct location p times. Because
#     w is sorted by (type_obs, xcoord, ycoord, j), a point's p rows are
#     consecutive, so the distinct locations are w[seq(1, nrow(w), by = p)] and
#     the full mask is rep(mask, each = p). The block structure is verified
#     before it is relied on (a cheap vectorised comparison, far cheaper than
#     the geometry it avoids), with a fallback to the original all-rows test if
#     it does not hold.
#
# (2) Xis[[i]] is by construction X[mark.pp[i] == X$marks], so its mask is just
#     the X mask subset by type. The p extra inside.owin() calls are replaced
#     by indexing.
#
# Total geometry work drops from n * (p + 2) to 2n tests -- a (p + 2) / 2 fold
# reduction, so 5x at p = 8. This matters for polygonal windows, where
# inside.owin() costs O(number of edges) per point; for rectangular windows the
# whole block is already negligible.
#
# `erodedwindow` may be supplied to skip the erosion() call entirely. This is
# worth doing in the grid search in SemiMarkov(), which calls
# SemiMarkov_fixed_R() once per (R_within, R_between, sat) combination and
# recomputes the identical erosion every time -- erosion() of a finely
# discretised polygon is not free (27 ms for a 4096-sided window here).
apply_edge_correction <- function(
  X, Xis, nis, w, p, mark.pp, edgecorrection, erodedwindow = NULL
) {
  if (is.null(erodedwindow)) {
    erodedwindow <- erosion(X$window, edgecorrection)
  }

  # --- (1) mask for w, testing each distinct location once ---
  nw <- nrow(w)
  fast_ok <- FALSE
  if (nw > 0L && nw %% p == 0L) {
    idx <- seq.int(1L, nw, by = p)
    ux <- w$xcoord[idx]
    uy <- w$ycoord[idx]
    fast_ok <- all(w$xcoord == rep(ux, each = p)) &&
      all(w$ycoord == rep(uy, each = p))
  }

  if (fast_ok) {
    keep_u <- inside.owin(x = ux, y = uy, w = erodedwindow)
    pts_in_window <- rep(keep_u, each = p)
  } else {
    pts_in_window <- inside.owin(
      x = w$xcoord,
      y = w$ycoord,
      w = erodedwindow
    )
  }
  w <- w[pts_in_window]

  # --- (2) one mask for X, reused for every Xis[[i]] ---
  keep_X <- inside.owin(X, w = erodedwindow)
  marks_int <- as.integer(X$marks)
  for (i in 1:p) {
    keep_i <- keep_X[marks_int == as.integer(mark.pp[i])]
    Xis[[i]] <- Xis[[i]][keep_i, ]
    nis[i] <- sum(keep_i)
  }
  X <- X[keep_X, ]

  list(
    X = X,
    Xis = Xis,
    nis = nis,
    w = w,
    erodedwindow = erodedwindow
  )
}
