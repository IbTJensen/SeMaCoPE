# Fast reimplementation of beta_estimation() and type_pred() -------------------
# Kept in its own file (rather than replacing the originals in SemiMarkov.R)
# so the two versions can be benchmarked and validated against each other,
# e.g.
#   identical(beta_estimation(w, betastart, p, q, Poisson)$par,
#             beta_estimation_fast(w, betastart, p, q, Poisson)$par)
#
# Both functions below avoid the merge.data.table() + order() calls that the
# originals perform on every optim() function evaluation. The grouping
# (xcoord, ycoord, type_obs) and the "true type" row indices do not change
# across evaluations, so they are computed once and reused, and the
# per-evaluation group sums use rowsum() instead of an aggregate-then-join.

# Fast replacement for type_pred(). Produces the same columns, column order
# and row order as the original, but replaces the aggregate + merge.data.table
# join with an in-place grouped assignment, and order() with setorder() to
# avoid the extra copy order() makes.
type_pred_fast <- function(w, beta) {
  cond_intensity <- as.vector(exp(as.matrix(w[, -(1:4), with = FALSE]) %*% beta))

  pts_lam <- w[, 1:4]
  pts_lam[, lambda := cond_intensity]
  pts_lam[, Lambda := sum(lambda), by = .(xcoord, ycoord, type_obs)]
  pts_lam[, prob := lambda / Lambda]
  data.table::setorder(pts_lam, type_obs, xcoord, ycoord, j)

  return(pts_lam)
}

# Fast replacement for beta_estimation(). Same signature, same optimizer
# (BFGS with a finite-difference gradient, same as the original) and the
# same betahat/convergence behaviour -- only the per-evaluation cost of
# loglik() changes.
beta_estimation_fast <- function(w, betastart, p, q, Poisson) {
  Xmat <- as.matrix(w[, -(1:4), with = FALSE])

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
    cond_intensity <- as.vector(exp(Xmat %*% beta))
    Lambda <- rowsum(cond_intensity, grp)[, 1]
    prob_true <- cond_intensity[true_idx] / Lambda[true_grp]
    sum(log(prob_true[!is.na(prob_true)]))
  }

  obj <- optim(
    par = betastart,
    fn = loglik,
    method = "BFGS",
    control = list(fnscale = -1)
  )
  return(obj)
}
