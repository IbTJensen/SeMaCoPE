// Compiled Newton-Raphson solver for beta_estimation() ------------------------
// Companion to R/beta_estimation_rcpp.R. Kept outside src/ (which R CMD
// build/install auto-compiles as part of the package DLL) so this file is
// purely opt-in: it is only touched when beta_estimation_rcpp() calls
// Rcpp::sourceCpp() on it, and never affects a normal package build/install.
//
// Computes loglik, gradient, and Sensitivity (= -Hessian, the same quantity
// already used elsewhere in this package as the Fisher-information-like "S"
// matrix, e.g. in Standard_error_matrix()) in a single fused pass over the
// design matrix per Newton iteration -- avoiding the redundant exp()/rowsum()
// recomputation that separate loglik()/score()/sensitivity() R closures each
// pay for independently. Softmax probabilities are computed with the
// standard per-group max-subtraction (log-sum-exp) trick for numerical
// stability, and each Newton step is halved until it doesn't decrease the
// log-likelihood -- both added because a raw, undamped Newton step from a
// naive implementation was observed to overflow exp() on the 2nd iteration
// (see beta_estimation()'s maxNR attempt).
//
// MISSING DATA: this routine requires a design matrix that is entirely
// finite. Covariates legitimately missing at some spatial locations are
// handled by beta_estimation_rcpp(), which drops whole groups (observed
// points) containing any non-finite covariate before calling in -- matching
// the complete-case-by-group semantics of the R beta_estimation(), where a
// single NA makes that group's Lambda (and hence every prob in the group) NA
// and the group drops out of `sum(log(prob[!is.na(prob)]))`. Passing a
// non-finite design matrix here is a hard error rather than silently
// producing a NaN-filled Hessian.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

namespace {

struct LGH {
  double loglik;
  arma::vec grad;
  arma::mat S;
};

// One fused pass: eta -> (stabilized) softmax probabilities -> loglik,
// gradient, and Sensitivity, all from the same intermediate quantities.
LGH loglik_grad_hess(
    const arma::mat& Wmat,
    const arma::ivec& grp0,       // 0-indexed group id per row
    const arma::uvec& true_idx0,  // 0-indexed row index of the true type, one per group, ordered by group
    const arma::vec& beta,
    int G
) {
  const int n = Wmat.n_rows;
  const int k = Wmat.n_cols;

  arma::vec eta = Wmat * beta;

  arma::vec group_max(G);
  group_max.fill(-arma::datum::inf);
  for (int i = 0; i < n; i++) {
    int g = grp0[i];
    if (eta[i] > group_max[g]) group_max[g] = eta[i];
  }

  arma::vec exp_shifted(n);
  arma::vec Lambda(G, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    int g = grp0[i];
    double v = std::exp(eta[i] - group_max[g]);
    exp_shifted[i] = v;
    Lambda[g] += v;
  }

  arma::vec prob(n);
  for (int i = 0; i < n; i++) {
    prob[i] = exp_shifted[i] / Lambda[grp0[i]];
  }

  // Per-group weighted mean of the design rows, E_prob_g[x].
  arma::mat Ex(G, k, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    Ex.row(grp0[i]) += prob[i] * Wmat.row(i);
  }

  double loglik = 0.0;
  for (int g = 0; g < G; g++) {
    double eta_true_g = eta[true_idx0[g]];
    loglik += eta_true_g - group_max[g] - std::log(Lambda[g]);
  }

  arma::mat Wtrue = Wmat.rows(true_idx0); // G x k
  arma::vec grad = (arma::sum(Wtrue, 0) - arma::sum(Ex, 0)).t();

  // S = sum_g (E_prob_g[x x'] - E_prob_g[x] E_prob_g[x]')
  //   = Wmat' diag(prob) Wmat - Ex' Ex
  arma::mat S = Wmat.t() * (Wmat.each_col() % prob) - Ex.t() * Ex;

  LGH out;
  out.loglik = loglik;
  out.grad = grad;
  out.S = S;
  return out;
}

} // namespace

// [[Rcpp::export]]
List beta_newton_cpp(
    const arma::mat& Wmat,
    const arma::ivec& grp,       // 1-indexed group id per row, contiguous 1..G
    const arma::ivec& true_idx,  // 1-indexed row of the true type, one per group, ordered by group
    arma::vec beta,
    int maxit = 100,
    double tol = 1e-8,
    int max_halving = 30,
    double ridge = 1e-8          // RELATIVE to mean(diag(S)), not absolute
) {
  const int n = Wmat.n_rows;
  const int k = Wmat.n_cols;

  // ---- Validation -----------------------------------------------------
  // These were previously assumed. A mismatch used to index out of bounds
  // (arma's operator[] is unchecked) and silently corrupt the Hessian, which
  // then surfaced only as an inscrutable "system is singular" warning.
  if ((int) grp.n_elem != n) {
    stop("beta_newton_cpp(): length(grp) must equal nrow(Wmat).");
  }
  if ((int) beta.n_elem != k) {
    stop("beta_newton_cpp(): length(beta) must equal ncol(Wmat).");
  }
  if (n == 0) {
    stop("beta_newton_cpp(): design matrix has no rows.");
  }
  if (!Wmat.is_finite()) {
    stop(
      "beta_newton_cpp(): design matrix contains non-finite values. Groups "
      "with missing covariates must be removed before calling; "
      "beta_estimation_rcpp() does this."
    );
  }
  if (grp.min() < 1) {
    stop("beta_newton_cpp(): group ids must be 1-indexed.");
  }

  const int G = grp.max();
  arma::ivec grp_count(G, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    grp_count[grp[i] - 1]++;
  }
  if (grp_count.min() == 0) {
    stop(
      "beta_newton_cpp(): group ids must be contiguous 1..G with at least "
      "one row per group."
    );
  }
  if ((int) true_idx.n_elem != G) {
    stop(
      "beta_newton_cpp(): expected exactly one 'true type' row per group, "
      "got " + std::to_string(true_idx.n_elem) + " for " +
      std::to_string(G) + " groups."
    );
  }
  if (true_idx.min() < 1 || true_idx.max() > n) {
    stop("beta_newton_cpp(): true_idx contains out-of-range row indices.");
  }

  arma::ivec grp0 = grp - 1;
  arma::uvec true_idx0 = arma::conv_to<arma::uvec>::from(true_idx - 1);
  for (int g = 0; g < G; g++) {
    if (grp0[true_idx0[g]] != g) {
      stop(
        "beta_newton_cpp(): true_idx must be ordered by group, with entry g "
        "pointing at a row belonging to group g."
      );
    }
  }

  // ---- Newton iteration ------------------------------------------------
  LGH cur = loglik_grad_hess(Wmat, grp0, true_idx0, beta, G);
  if (!std::isfinite(cur.loglik)) {
    stop("beta_newton_cpp(): log-likelihood is non-finite at the start value.");
  }

  int iter = 0;
  bool converged = false;
  bool step_failed = false;
  std::string failure_msg;

  for (iter = 1; iter <= maxit; iter++) {
    // Ridge is scaled by the magnitude of S: entries of S grow with the
    // number of points, so a fixed absolute ridge is numerically invisible
    // on real-sized problems and does nothing to regularise a rank-deficient
    // system.
    double s_scale = arma::trace(cur.S) / (double) k;
    if (!std::isfinite(s_scale) || s_scale <= 0.0) s_scale = 1.0;

    arma::vec step;
    bool solved = false;
    // Escalating ridge. solve_opts::no_approx makes a singular system return
    // false instead of emitting arma's "attempting approx solution" warning
    // and handing back an unusable step, so the escalation is under our
    // control rather than arma's.
    const double mult[4] = {0.0, 1.0, 1e2, 1e4};
    for (int t = 0; t < 4 && !solved; t++) {
      arma::mat Sreg = cur.S;
      if (t > 0) Sreg.diag() += ridge * mult[t] * s_scale;
      solved = arma::solve(step, Sreg, cur.grad, arma::solve_opts::no_approx);
      if (solved && !step.is_finite()) solved = false;
    }
    if (!solved) {
      step_failed = true;
      failure_msg = "Hessian could not be factorised even with ridge "
                    "regularisation (design may be rank-deficient).";
      break;
    }

    double step_size = 1.0;
    arma::vec beta_new;
    LGH next;
    bool accepted = false;
    for (int h = 0; h < max_halving; h++) {
      beta_new = beta + step_size * step;
      next = loglik_grad_hess(Wmat, grp0, true_idx0, beta_new, G);
      if (std::isfinite(next.loglik) && next.loglik >= cur.loglik - 1e-10) {
        accepted = true;
        break;
      }
      step_size *= 0.5;
    }
    if (!accepted) {
      step_failed = true;
      failure_msg = "no step size improved the log-likelihood after " +
                    std::to_string(max_halving) + " halvings.";
      break;
    }

    double step_norm = arma::norm(step_size * step, 2);
    beta = beta_new;
    cur = next;

    if (step_norm < tol) {
      converged = true;
      break;
    }
  }

  if (!converged && !step_failed) {
    failure_msg = "reached maxit (" + std::to_string(maxit) +
                  ") without meeting the convergence tolerance.";
  }

  return List::create(
    Named("par") = beta,
    Named("value") = cur.loglik,
    Named("gradient") = cur.grad,
    Named("Sensitivity") = cur.S,
    Named("iterations") = iter,
    Named("convergence") = converged ? 0 : 1,
    Named("message") = converged ? "" : failure_msg
  );
}
