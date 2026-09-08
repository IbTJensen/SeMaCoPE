// Compiled core of Covariate_setup(), with separate within/between saturation -
// Companion to R/Covariate_setup_sat2_rcpp.R.
//
// WHAT DIFFERS FROM covariate_setup_core_cpp()
// --------------------------------------------
// Nothing except where the saturation parameter comes from. The existing core
// takes `sat_i`, one value per *type* l (in practice always rep(sat, p), a
// single global saturation). This version takes two scalars and applies them
// exactly the way the two interaction radii are already applied:
//
//   R  [m,l] = R_within   if m == l, else R_between      (unchanged)
//   sat[m,l] = sat_within if m == l, else sat_between    (new)
//
// so a within-type interaction saturates at sat_within and a between-type one
// at sat_between. Setting sat_within == sat_between reproduces
// covariate_setup_core_cpp() with sat_i = rep(sat, p) exactly. Either value
// may be Inf independently, so e.g. sat_within = Inf with a finite
// sat_between gives Strauss within types and Geyer saturation between them.
//
// The quantities computed are therefore
//
//   s      _{k,l}(u) = min( N_{k,l}(u), sat[k,l] )
//   DeltaS_{m,l}(u) = # { v : type(v) = m, d(u,v) <= R[m,l],
//                             N_{m,l}(v) - [type(u) == l] <= sat[m,l] - 1 }
//   N_{m,l}(v)      = # { z != v : type(z) = l, d(v,z) <= R[m,l] }
//
// The `- [type(u) == l]` term is the leave-u-out correction: adding a type-l
// point at u raises neighbour v's saturated type-l count only if v's count on
// the pattern *with u removed* is still below sat, and N is computed over the
// full observed pattern. See src/covariate_setup_rcpp.cpp for the full
// discussion; the logic is identical here.
//
// NOTE ON DUPLICATION
// -------------------
// The cell-grid helpers below are a copy of those in
// src/covariate_setup_rcpp.cpp rather than a shared header, so that adding
// this variant leaves the existing, validated translation unit untouched.
// They live in an anonymous namespace (internal linkage), so the two copies
// cannot collide at link time. Worth consolidating into a shared header if
// this variant is kept.

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

namespace {

// Uniform cell list over a subset of points, cell size = query radius, so a
// radius query never has to look beyond the 3x3 block of cells around the
// query point.
struct CellGrid {
  double cell = 0.0, xmin = 0.0, ymin = 0.0;
  int nx = 0, ny = 0;
  std::vector<int> start;  // CSR offsets, length nx*ny + 1
  std::vector<int> items;  // global point indices
  bool empty = true;
};

void grid_build(
    CellGrid& g,
    const std::vector<int>& idx,
    const double* x,
    const double* y,
    double radius
) {
  g.empty = idx.empty() || !(radius > 0.0) || !R_finite(radius);
  if (g.empty) return;

  g.cell = radius;
  double xmn = x[idx[0]], xmx = xmn, ymn = y[idx[0]], ymx = ymn;
  for (std::size_t t = 1; t < idx.size(); ++t) {
    const int id = idx[t];
    if (x[id] < xmn) xmn = x[id];
    if (x[id] > xmx) xmx = x[id];
    if (y[id] < ymn) ymn = y[id];
    if (y[id] > ymx) ymx = y[id];
  }
  g.xmin = xmn;
  g.ymin = ymn;
  g.nx = (int) std::floor((xmx - xmn) / g.cell) + 1;
  g.ny = (int) std::floor((ymx - ymn) / g.cell) + 1;
  if (g.nx < 1) g.nx = 1;
  if (g.ny < 1) g.ny = 1;

  const std::size_t ncell = (std::size_t) g.nx * (std::size_t) g.ny;
  g.start.assign(ncell + 1, 0);
  std::vector<int> cid(idx.size());
  for (std::size_t t = 0; t < idx.size(); ++t) {
    const int id = idx[t];
    int cx = (int) std::floor((x[id] - g.xmin) / g.cell);
    int cy = (int) std::floor((y[id] - g.ymin) / g.cell);
    if (cx < 0) cx = 0;
    if (cx >= g.nx) cx = g.nx - 1;
    if (cy < 0) cy = 0;
    if (cy >= g.ny) cy = g.ny - 1;
    const int c = cy * g.nx + cx;
    cid[t] = c;
    g.start[c + 1]++;
  }
  for (std::size_t c = 1; c <= ncell; ++c) g.start[c] += g.start[c - 1];

  g.items.assign(idx.size(), 0);
  std::vector<int> pos(g.start.begin(), g.start.begin() + ncell);
  for (std::size_t t = 0; t < idx.size(); ++t) g.items[pos[cid[t]]++] = idx[t];
}

// Sum of wts[v] over grid points v with d((qx,qy), v) <= radius, skipping
// v == self_idx. wts == nullptr means unit weights (a plain count).
double grid_query(
    const CellGrid& g,
    double qx, double qy, double radius,
    const double* x, const double* y,
    const double* wts, int self_idx
) {
  if (g.empty) return 0.0;
  const double r2 = radius * radius;
  const int cx = (int) std::floor((qx - g.xmin) / g.cell);
  const int cy = (int) std::floor((qy - g.ymin) / g.cell);
  double acc = 0.0;
  for (int dy = -1; dy <= 1; ++dy) {
    const int yy = cy + dy;
    if (yy < 0 || yy >= g.ny) continue;
    for (int dx = -1; dx <= 1; ++dx) {
      const int xx = cx + dx;
      if (xx < 0 || xx >= g.nx) continue;
      const int c = yy * g.nx + xx;
      for (int t = g.start[c]; t < g.start[c + 1]; ++t) {
        const int id = g.items[t];
        if (id == self_idx) continue;
        const double ddx = x[id] - qx, ddy = y[id] - qy;
        if (ddx * ddx + ddy * ddy <= r2) acc += (wts ? wts[id] : 1.0);
      }
    }
  }
  return acc;
}

} // namespace

// type: 1-indexed type code per point (1..p).
// sat_within / sat_between: saturation for within-type and between-type
//   interactions respectively; either may be Inf.
// Returns s and DeltaS as n x p^2 matrices and Int_mat as (n*p) x p^2, all
// with columns in the (k outer, l inner) order the R code assumes, i.e.
// column (k-1)*p + l holds the "k-l" term.
// [[Rcpp::export]]
List covariate_setup_core_sat2_cpp(
    NumericVector x,
    NumericVector y,
    IntegerVector type,
    int p,
    double R_within,
    double R_between,
    double sat_within,
    double sat_between,
    bool want_components = true
) {
  const int n = x.size();
  if (y.size() != n || type.size() != n) {
    stop("covariate_setup_core_sat2_cpp(): x, y and type must have equal length.");
  }
  if (p < 1) stop("covariate_setup_core_sat2_cpp(): p must be >= 1.");
  if (!R_finite(R_within) || !R_finite(R_between) ||
      R_within < 0 || R_between < 0) {
    stop("covariate_setup_core_sat2_cpp(): R_within and R_between must be finite and non-negative.");
  }
  // Inf is allowed (it is what makes an interaction Strauss rather than
  // Geyer), NaN and negatives are not.
  if (ISNAN(sat_within) || ISNAN(sat_between) ||
      sat_within < 0 || sat_between < 0) {
    stop("covariate_setup_core_sat2_cpp(): sat_within and sat_between must be non-negative (Inf allowed).");
  }
  for (int i = 0; i < n; ++i) {
    if (type[i] < 1 || type[i] > p) {
      stop("covariate_setup_core_sat2_cpp(): type codes must lie in 1..p.");
    }
  }

  const double* px = x.begin();
  const double* py = y.begin();
  const int p2 = p * p;

  std::vector< std::vector<int> > idx_by_type(p);
  for (int i = 0; i < n; ++i) idx_by_type[type[i] - 1].push_back(i);

  // R[m,l] takes only two values, so 2p grids cover every query in both
  // passes below.
  std::vector<CellGrid> grid_w(p), grid_b(p);
  for (int l = 0; l < p; ++l) {
    grid_build(grid_w[l], idx_by_type[l], px, py, R_within);
    grid_build(grid_b[l], idx_by_type[l], px, py, R_between);
  }

  // Pass 1: N_within[l][v] / N_between[l][v] = number of type-(l+1) points
  // within the respective radius of v, excluding v itself. Unaffected by the
  // saturation parameters.
  std::vector< std::vector<double> > Nw(p, std::vector<double>(n, 0.0));
  std::vector< std::vector<double> > Nb(p, std::vector<double>(n, 0.0));
  for (int l = 0; l < p; ++l) {
    for (int i = 0; i < n; ++i) {
      Nw[l][i] = grid_query(grid_w[l], px[i], py[i], R_within, px, py, nullptr, i);
      Nb[l][i] = grid_query(grid_b[l], px[i], py[i], R_between, px, py, nullptr, i);
    }
    Rcpp::checkUserInterrupt();
  }

  // s_{k,l}(v) = min(N_{k,l}(v), sat[k,l]), with both the count radius and
  // the saturation chosen by whether the interaction is within or between.
  NumericMatrix s(n, p2);
  for (int k = 0; k < p; ++k) {
    for (int l = 0; l < p; ++l) {
      const bool same = (k == l);
      const std::vector<double>& N = same ? Nw[l] : Nb[l];
      const double cap = same ? sat_within : sat_between;
      const int c = k * p + l;
      for (int i = 0; i < n; ++i) s(i, c) = (N[i] > cap) ? cap : N[i];
    }
  }

  // Pass 2: DeltaS_{m,l}(u). Two weight vectors selected by the query point's
  // own type -- the saturation indicator is a property of the pair (u, v),
  // not of v alone.
  NumericMatrix DeltaS(n, p2);
  std::vector<double> wt_other(n, 0.0), wt_self(n, 0.0);
  for (int m = 0; m < p; ++m) {
    for (int l = 0; l < p; ++l) {
      const bool same = (m == l);
      const std::vector<double>& N = same ? Nw[l] : Nb[l];
      // Inf - 1 == Inf, so a saturation of Inf keeps every neighbour.
      const double thresh = (same ? sat_within : sat_between) - 1.0;
      const double radius = same ? R_within : R_between;
      const CellGrid& g = same ? grid_w[m] : grid_b[m];

      std::fill(wt_other.begin(), wt_other.end(), 0.0);
      std::fill(wt_self.begin(), wt_self.end(), 0.0);
      for (std::size_t t = 0; t < idx_by_type[m].size(); ++t) {
        const int v = idx_by_type[m][t];
        wt_other[v] = (N[v] <= thresh) ? 1.0 : 0.0;
        wt_self[v] = (N[v] - 1.0 <= thresh) ? 1.0 : 0.0;
      }

      const int c = m * p + l;
      for (int u = 0; u < n; ++u) {
        const double* wsel =
          (type[u] - 1 == l) ? wt_self.data() : wt_other.data();
        DeltaS(u, c) = grid_query(g, px[u], py[u], radius, px, py, wsel, u);
      }
      Rcpp::checkUserInterrupt();
    }
  }

  // Int_mat: rows stacked by candidate type j (all n points for j = 1, then
  // j = 2, ...), matching Reduce(rbind, Int_mat_list) and the `pts` table.
  // For candidate j, the "k-l" column takes s_{j,l} when k == j and
  // DeltaS_{k,j} when l == j (the diagonal "j-j" column takes both).
  NumericMatrix Int_mat(n * p, p2);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) {
      const int row = j * n + i;
      for (int k = 0; k < p; ++k) {
        for (int l = 0; l < p; ++l) {
          const int c = k * p + l;
          double val = 0.0;
          if (k == j) val += s(i, c);
          if (l == j) val += DeltaS(i, c);
          Int_mat(row, c) = val;
        }
      }
    }
  }

  if (!want_components) {
    return List::create(Named("Int_mat") = Int_mat);
  }
  return List::create(
    Named("Int_mat") = Int_mat,
    Named("s") = s,
    Named("DeltaS") = DeltaS
  );
}
