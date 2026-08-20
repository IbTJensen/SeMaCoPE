#' SeMaCoPE: Semi-parametric Markov models for multi-type point patterns
#'
#' Implementation of a semi-parametric Markov model for multi-type point
#' pattern data, fitted by conditional pseudo-likelihood.
#'
#' The computational core is in C++ (see `src/`): the design-matrix
#' construction in `covariate_setup_core_cpp()` and the Newton-Raphson solver
#' in `beta_newton_cpp()` / `beta_newton_omp_cpp()`. These are compiled into
#' the package library at install time and reached through the generated
#' `RcppExports` wrappers, so nothing depends on the working directory.
#'
#' @useDynLib SeMaCoPE, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom data.table data.table as.data.table rbindlist setcolorder setorder melt dcast merge.data.table
#' @importFrom spatstat.geom closepairs crosspairs erosion inside.owin is.im ppp
#' @importFrom Matrix sparseMatrix
#' @importFrom VGAM vglm multinomial
#' @importFrom stats optim model.matrix predict sd pnorm na.omit runif coef
#' @importFrom utils globalVariables
#' @keywords internal
"_PACKAGE"

# data.table's NSE (`:=`, `.N`, `.GRP`, `.SD`, and bare column names) is
# invisible to R CMD check's global-variable analysis. Declaring them here
# keeps the check output clean without changing behaviour.
globalVariables(c(
  ".", ".N", ".GRP", ".SD", ":=",
  "xcoord", "ycoord", "type_obs", "j", "l", "l2", "grp",
  "lambda", "Lambda", "prob", "sat", "a",
  "Neighbour_x", "Neighbour_y", "Neighbour_type",
  "s_Strauss_kl_v", "DeltaS_kl", "kl_name",
  # data.table's `..name` (look up `name` in the calling frame) syntax
  "..between_idx", "..covar_cols", "..Covariate_cols", "..Intercept_cols",
  "..m", "..p", "..rep_col", "..var_cols", "..within_idx"
))
