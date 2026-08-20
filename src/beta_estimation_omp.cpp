// Multi-threaded (OpenMP) Newton-Raphson solver for beta_estimation() ---------
// Companion to R/beta_estimation_omp.R; parallel counterpart of
// src/beta_estimation_rcpp.cpp. Compiled into the package DLL by
// R CMD INSTALL; see src/Makevars for the OpenMP flags.
//
// WHY THIS IS WORTH PARALLELISING
// -------------------------------
// Profiling the serial kernel shows the Sensitivity matrix
//   S = Wmat' diag(prob) Wmat - Ex' Ex
// is 70-76% of per-iteration time and scales as O(n k^2), so it dominates
// exactly in the regime that motivated all of this (large numbers of
// covariates). Measured on the serial version, 6 candidate types, 2000
// points: 16ms/iteration at k=40 rising to 1798ms/iteration at k=400.
//
// R here links against the reference BLAS, which is single-threaded, so that
// dgemm gets no parallelism for free. (With a threaded BLAS -- OpenBLAS, MKL
// -- it partly would, and stacking OpenMP on top risks oversubscription; see
// the `nthreads` note in the R wrapper.)
//
// TWO INDEPENDENT WINS
// --------------------
// 1. Symmetry. prob is a softmax probability, hence >= 0, so
//      Wmat' diag(prob) Wmat = B'B  with  B = sqrt(prob) .* Wmat,
//    which is a symmetric rank-k update: only the upper triangle needs
//    computing, roughly halving the flops before any threading.
// 2. Threading. The upper triangle is computed in column blocks, each block
//    writing a disjoint set of columns of S, so no locks or per-thread k x k
//    accumulators are needed. Everything else (group max, Lambda, prob, Ex,
//    per-group log-likelihood) is parallelised over *groups* using a CSR
//    index, so each thread touches only its own groups' rows and its own row
//    of Ex -- again conflict-free.
//
// DETERMINISM
// -----------
// Results are reproducible: no floating-point reduction depends on thread
// scheduling. Per-group log-likelihood contributions are written to an array
// and summed serially in fixed order rather than via `reduction(+:)`, and all
// matrix writes are to disjoint memory. Output is therefore identical run to
// run and independent of nthreads. It is NOT bit-identical to
// beta_newton_cpp(), because the B'B formulation rounds differently from
// Wmat'(Wmat .* prob); agreement is to ~1e-12 relative.

// [[Rcpp::depends(RcppArmadillo)]]
// NOTE: no `// [[Rcpp::plugins(openmp)]]` here -- that attribute is honoured
// only by Rcpp::sourceCpp(). Inside a package the OpenMP flags come from
// src/Makevars; the #ifdef _OPENMP guards below keep this file compiling (and
// running, serially) on toolchains without OpenMP.
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

namespace {

struct LGH {
  double loglik;
  arma::vec grad;
  arma::mat S;
};

// Reused across Newton iterations so the per-iteration cost is arithmetic
// only, with no repeated allocation of the n x k temporaries.
struct Workspace {
  arma::vec eta, expv, prob, gmax, lam, ll_g;
  arma::mat Ex, B;
  void init(int n, int k, int G) {
    eta.set_size(n); expv.set_size(n); prob.set_size(n);
    gmax.set_size(G); lam.set_size(G); ll_g.set_size(G);
    Ex.set_size(G, k); B.set_size(n, k);
  }
};

// C = A'A (k x k, symmetric PSD), upper triangle computed in column blocks in
// parallel, then mirrored. Column ranges of A are contiguous in column-major
// storage, so lightweight non-owning views avoid copying any submatrix.
void syrk_upper_parallel(const arma::mat& A, arma::mat& C, int nt, int block) {
  const int k = (int) A.n_cols;
  const int n = (int) A.n_rows;
  C.zeros(k, k);
  const int nblk = (k + block - 1) / block;
  double* Aptr = const_cast<double*>(A.memptr());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nt)
#endif
  for (int b = 0; b < nblk; ++b) {
    const int j0 = b * block;
    const int j1 = std::min(k, j0 + block);
    const arma::mat left(Aptr, n, j1, false, true);              // cols 0..j1-1
    const arma::mat right(Aptr + (std::size_t) j0 * n, n, j1 - j0, false, true);
    C.submat(0, j0, j1 - 1, j1 - 1) = left.t() * right;
  }

  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < i; ++j) C(i, j) = C(j, i);
  }
}

LGH loglik_grad_hess_omp(
    const arma::mat& Wmat,
    const std::vector<int>& g_start,   // CSR offsets over groups, length G+1
    const std::vector<int>& g_rows,    // row indices grouped by group
    const arma::uvec& true_idx0,
    const arma::vec& beta,
    int G, int nt, int block,
    Workspace& ws
) {
  const int n = (int) Wmat.n_rows;
  const int k = (int) Wmat.n_cols;

  ws.eta = Wmat * beta;

  const double* etap = ws.eta.memptr();
  double* expvp = ws.expv.memptr();
  double* probp = ws.prob.memptr();
  double* llg = ws.ll_g.memptr();
  ws.Ex.zeros();

  // Parallel over groups: group g touches only its own rows and Ex.row(g).
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (int g = 0; g < G; ++g) {
    const int b0 = g_start[g], b1 = g_start[g + 1];

    double gmax = -arma::datum::inf;
    for (int t = b0; t < b1; ++t) {
      const double e = etap[g_rows[t]];
      if (e > gmax) gmax = e;
    }

    double lam = 0.0;
    for (int t = b0; t < b1; ++t) {
      const int r = g_rows[t];
      const double v = std::exp(etap[r] - gmax);
      expvp[r] = v;
      lam += v;
    }

    double* exrow = ws.Ex.colptr(0) + g;  // Ex is G x k, column-major: stride G
    for (int t = b0; t < b1; ++t) {
      const int r = g_rows[t];
      const double pr = expvp[r] / lam;
      probp[r] = pr;
      for (int c = 0; c < k; ++c) exrow[(std::size_t) c * G] += pr * Wmat(r, c);
    }

    llg[g] = etap[true_idx0[g]] - gmax - std::log(lam);
  }

  // Summed serially in fixed order so the result does not depend on the
  // thread schedule.
  double loglik = 0.0;
  for (int g = 0; g < G; ++g) loglik += llg[g];

  // B = sqrt(prob) .* Wmat, so that Wmat' diag(prob) Wmat == B'B (syrk).
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (int c = 0; c < k; ++c) {
    double* bcol = ws.B.colptr(c);
    const double* wcol = Wmat.colptr(c);
    for (int r = 0; r < n; ++r) bcol[r] = std::sqrt(probp[r]) * wcol[r];
  }

  arma::mat S;
  syrk_upper_parallel(ws.B, S, nt, block);

  arma::mat ExtEx;
  syrk_upper_parallel(ws.Ex, ExtEx, nt, block);
  S -= ExtEx;

  arma::mat Wtrue = Wmat.rows(true_idx0);
  arma::vec grad = (arma::sum(Wtrue, 0) - arma::sum(ws.Ex, 0)).t();

  LGH out;
  out.loglik = loglik;
  out.grad = grad;
  out.S = S;
  return out;
}

} // namespace

// [[Rcpp::export]]
int omp_max_threads_cpp() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

// [[Rcpp::export]]
List beta_newton_omp_cpp(
    const arma::mat& Wmat,
    const arma::ivec& grp,       // 1-indexed group id per row, contiguous 1..G
    const arma::ivec& true_idx,  // 1-indexed row of the true type, one per group, ordered by group
    arma::vec beta,
    int nthreads = 1,
    int maxit = 100,
    double tol = 1e-8,
    int max_halving = 30,
    double ridge = 1e-8,         // RELATIVE to mean(diag(S)), not absolute
    int block = 64
) {
  const int n = (int) Wmat.n_rows;
  const int k = (int) Wmat.n_cols;

  // ---- Validation (identical to beta_newton_cpp) -----------------------
  if ((int) grp.n_elem != n) stop("beta_newton_omp_cpp(): length(grp) must equal nrow(Wmat).");
  if ((int) beta.n_elem != k) stop("beta_newton_omp_cpp(): length(beta) must equal ncol(Wmat).");
  if (n == 0) stop("beta_newton_omp_cpp(): design matrix has no rows.");
  if (!Wmat.is_finite()) {
    stop(
      "beta_newton_omp_cpp(): design matrix contains non-finite values. Groups "
      "with missing covariates must be removed before calling; "
      "beta_estimation_omp() does this."
    );
  }
  if (grp.min() < 1) stop("beta_newton_omp_cpp(): group ids must be 1-indexed.");

  const int G = grp.max();
  arma::ivec grp_count(G, arma::fill::zeros);
  for (int i = 0; i < n; i++) grp_count[grp[i] - 1]++;
  if (grp_count.min() == 0) {
    stop("beta_newton_omp_cpp(): group ids must be contiguous 1..G with at least one row per group.");
  }
  if ((int) true_idx.n_elem != G) {
    stop(
      "beta_newton_omp_cpp(): expected exactly one 'true type' row per group, got " +
      std::to_string(true_idx.n_elem) + " for " + std::to_string(G) + " groups."
    );
  }
  if (true_idx.min() < 1 || true_idx.max() > n) {
    stop("beta_newton_omp_cpp(): true_idx contains out-of-range row indices.");
  }

  arma::ivec grp0 = grp - 1;
  arma::uvec true_idx0 = arma::conv_to<arma::uvec>::from(true_idx - 1);
  for (int g = 0; g < G; g++) {
    if (grp0[true_idx0[g]] != g) {
      stop("beta_newton_omp_cpp(): true_idx must be ordered by group, with entry g pointing at a row belonging to group g.");
    }
  }

  if (block < 1) block = 64;
  int nt = nthreads;
#ifdef _OPENMP
  if (nt < 1) nt = omp_get_max_threads();
  if (nt > omp_get_max_threads()) nt = omp_get_max_threads();
#else
  nt = 1;
#endif

  // CSR index of rows per group -- built once, lets the group loop be
  // parallelised without assuming rows of a group are contiguous.
  std::vector<int> g_start((std::size_t) G + 1, 0), g_rows((std::size_t) n, 0);
  for (int i = 0; i < n; ++i) g_start[grp0[i] + 1]++;
  for (int g = 0; g < G; ++g) g_start[g + 1] += g_start[g];
  {
    std::vector<int> pos(g_start.begin(), g_start.begin() + G);
    for (int i = 0; i < n; ++i) g_rows[pos[grp0[i]]++] = i;
  }

  Workspace ws;
  ws.init(n, k, G);

  LGH cur = loglik_grad_hess_omp(Wmat, g_start, g_rows, true_idx0, beta, G, nt, block, ws);
  if (!std::isfinite(cur.loglik)) {
    stop("beta_newton_omp_cpp(): log-likelihood is non-finite at the start value.");
  }

  int iter = 0;
  bool converged = false, step_failed = false;
  std::string failure_msg;

  for (iter = 1; iter <= maxit; iter++) {
    double s_scale = arma::trace(cur.S) / (double) k;
    if (!std::isfinite(s_scale) || s_scale <= 0.0) s_scale = 1.0;

    arma::vec step;
    bool solved = false;
    const double mult[4] = {0.0, 1.0, 1e2, 1e4};
    for (int t = 0; t < 4 && !solved; t++) {
      arma::mat Sreg = cur.S;
      if (t > 0) Sreg.diag() += ridge * mult[t] * s_scale;
      solved = arma::solve(step, Sreg, cur.grad, arma::solve_opts::no_approx);
      if (solved && !step.is_finite()) solved = false;
    }
    if (!solved) {
      step_failed = true;
      failure_msg = "Hessian could not be factorised even with ridge regularisation (design may be rank-deficient).";
      break;
    }

    double step_size = 1.0;
    arma::vec beta_new;
    LGH next;
    bool accepted = false;
    for (int h = 0; h < max_halving; h++) {
      beta_new = beta + step_size * step;
      next = loglik_grad_hess_omp(Wmat, g_start, g_rows, true_idx0, beta_new, G, nt, block, ws);
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
    if (step_norm < tol) { converged = true; break; }
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
    Named("message") = converged ? "" : failure_msg,
    Named("threads_used") = nt
  );
}
