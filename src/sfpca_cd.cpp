// Copyright (c) 2025 genpca contributors
//
// Coordinate-descent solver for the SFPCA rank-1 subproblem
//
//   min_x  0.5 * x' S x - b' x + P(x; lambda)
//
// where S = I + alpha * Omega is sparse SPD (column-compressed) and P is
// either the L1 penalty (lambda * ||x||_1) or the SCAD penalty. The solver
// maintains the gradient g = S x incrementally (each coordinate update costs
// one sparse column pass) and uses a glmnet-style active-set strategy:
// full sweeps alternate with sweeps restricted to the current support until
// a full sweep changes nothing.
#include <RcppEigen.h>
#include <cmath>
#include <vector>
// [[Rcpp::depends(RcppEigen)]]

namespace {

inline double soft_threshold(const double z, const double t) {
  if (z > t) return z - t;
  if (z < -t) return z + t;
  return 0.0;
}

// SCAD penalty value at |x| (includes the lambda factor).
inline double scad_value(const double ax, const double lambda, const double a) {
  if (ax <= lambda) return lambda * ax;
  if (ax <= a * lambda) {
    return (2.0 * a * lambda * ax - ax * ax - lambda * lambda) / (2.0 * (a - 1.0));
  }
  return (a + 1.0) * lambda * lambda / 2.0;
}

// Univariate SCAD-penalized quadratic:
//   min_x 0.5 * s * x^2 - z * x + scad(|x|; lambda, a),  s > 0.
// Solved by evaluating the stationary point of each SCAD region (clipped to
// its region) plus 0, and taking the minimizer. Robust at region boundaries
// and for any s > 0.
double scad_univariate(const double z, const double s, const double lambda,
                       const double a) {
  if (z == 0.0) return 0.0;
  const double sgn = (z > 0.0) ? 1.0 : -1.0;
  const double az = std::fabs(z);

  double best_x = 0.0;
  double best_h = 0.0;  // h(0) = 0

  auto consider = [&](double x) {
    if (x <= 0.0) return;
    const double h = 0.5 * s * x * x - az * x + scad_value(x, lambda, a);
    if (h < best_h) {
      best_h = h;
      best_x = x;
    }
  };

  // Region 1: 0 < x <= lambda, stationary point (az - lambda) / s
  {
    double x1 = (az - lambda) / s;
    if (x1 > 0.0) consider(std::min(x1, lambda));
  }
  // Region 2: lambda < x <= a*lambda, stationary point
  // (az - a*lambda/(a-1)) / (s - 1/(a-1)) when the curvature stays positive
  {
    const double denom = s - 1.0 / (a - 1.0);
    if (denom > 0.0) {
      double x2 = (az - a * lambda / (a - 1.0)) / denom;
      x2 = std::max(lambda, std::min(x2, a * lambda));
      consider(x2);
    } else {
      // Concave region: minimum at a boundary; both are checked by the
      // clipped candidates of regions 1 and 3, but check the far edge too.
      consider(a * lambda);
    }
  }
  // Region 3: x > a*lambda, stationary point az / s
  {
    double x3 = az / s;
    if (x3 > a * lambda) consider(x3);
  }

  return sgn * best_x;
}

}  // namespace

//' Coordinate descent for the SFPCA penalized quadratic subproblem
//'
//' Internal solver for `min_x 0.5 x'Sx - b'x + P(x; lambda)` with sparse SPD
//' `S` and an L1 or SCAD penalty.
//'
//' @param S sparse SPD matrix (`dgCMatrix`)
//' @param b numeric vector, linear term
//' @param x0 numeric vector, warm start
//' @param lambda penalty level (must be >= 0)
//' @param penalty 0 for L1, 1 for SCAD
//' @param scad_a SCAD shape parameter (> 2)
//' @param max_sweeps maximum number of full-equivalent sweeps
//' @param tol convergence tolerance on the KKT residual (gradient units)
//' @return list with `x`, `sweeps`, and `kkt` (max KKT residual)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List sfpca_cd_solve_cpp(const Eigen::Map<Eigen::SparseMatrix<double>> S,
                              const Eigen::Map<Eigen::VectorXd> b,
                              const Eigen::Map<Eigen::VectorXd> x0,
                              const double lambda, const int penalty,
                              const double scad_a, const int max_sweeps,
                              const double tol) {
  const int n = S.rows();
  if (S.cols() != n) Rcpp::stop("S must be square.");
  if (b.size() != n || x0.size() != n) Rcpp::stop("b and x0 must have length nrow(S).");
  if (!std::isfinite(lambda) || lambda < 0.0) Rcpp::stop("lambda must be non-negative and finite.");
  if (penalty != 0 && penalty != 1) Rcpp::stop("penalty must be 0 (l1) or 1 (scad).");
  if (penalty == 1 && scad_a <= 2.0) Rcpp::stop("scad_a must be > 2.");
  // NaN/Inf must be rejected up front: soft_threshold(NaN) returns 0 and
  // std::max(0.0, NaN) returns 0, so non-finite input would otherwise come
  // back as a "converged" solution with KKT residual 0.
  if (!b.allFinite()) Rcpp::stop("b must contain only finite values.");
  if (!x0.allFinite()) Rcpp::stop("x0 must contain only finite values.");
  for (int j = 0; j < S.outerSize(); ++j) {
    for (Eigen::Map<Eigen::SparseMatrix<double>>::InnerIterator it(S, j); it; ++it) {
      if (!std::isfinite(it.value())) Rcpp::stop("S must contain only finite values.");
    }
  }

  Eigen::VectorXd x = x0;
  Eigen::VectorXd diag(n);
  for (int j = 0; j < n; ++j) {
    double d = S.coeff(j, j);
    if (!std::isfinite(d) || d <= 0.0) Rcpp::stop("S must have strictly positive diagonal.");
    diag[j] = d;
  }

  // Gradient state g = S x
  Eigen::VectorXd g = S * x;

  // One coordinate update; returns S_jj * |delta| (gradient-scale change).
  auto update_coord = [&](const int j) -> double {
    const double xj_old = x[j];
    // z_j = b_j - sum_{k != j} S_jk x_k
    const double z = b[j] - (g[j] - diag[j] * xj_old);
    double xj_new;
    if (penalty == 0) {
      xj_new = soft_threshold(z, lambda) / diag[j];
    } else {
      xj_new = scad_univariate(z, diag[j], lambda, scad_a);
    }
    const double delta = xj_new - xj_old;
    if (delta != 0.0) {
      x[j] = xj_new;
      for (Eigen::Map<Eigen::SparseMatrix<double>>::InnerIterator it(S, j); it;
           ++it) {
        g[it.row()] += delta * it.value();
      }
    }
    return diag[j] * std::fabs(delta);
  };

  // KKT residual (subgradient optimality) from the maintained gradient g,
  // using the penalty derivative. Gradient units, so `tol` is compared
  // directly against it.
  auto kkt_residual = [&]() -> double {
    double kkt = 0.0;
    for (int j = 0; j < n; ++j) {
      const double grad = g[j] - b[j];
      double r;
      if (x[j] == 0.0) {
        // Threshold at zero is lambda for both L1 and SCAD.
        r = std::max(0.0, std::fabs(grad) - lambda);
      } else {
        const double ax = std::fabs(x[j]);
        double dp;  // penalty derivative at |x_j|
        if (penalty == 0) {
          dp = lambda;
        } else if (ax <= lambda) {
          dp = lambda;
        } else if (ax <= scad_a * lambda) {
          dp = (scad_a * lambda - ax) / (scad_a - 1.0);
        } else {
          dp = 0.0;
        }
        r = std::fabs(grad + ((x[j] > 0.0) ? dp : -dp));
      }
      kkt = std::max(kkt, r);
    }
    return kkt;
  };

  // Inner (active-set) sweeps stop once coordinate updates are small on the
  // gradient scale; the authoritative convergence test is the KKT residual
  // after each full sweep.
  const double inner_tol = 0.1 * tol;

  int sweeps = 0;
  bool converged = false;
  double kkt = 0.0;
  while (sweeps < max_sweeps) {
    Rcpp::checkUserInterrupt();
    // Full sweep
    for (int j = 0; j < n; ++j) update_coord(j);
    ++sweeps;
    kkt = kkt_residual();
    if (kkt < tol) {
      converged = true;
      break;
    }

    // Active-set sweeps over the current support
    std::vector<int> active;
    active.reserve(n);
    for (int j = 0; j < n; ++j) {
      if (x[j] != 0.0) active.push_back(j);
    }
    while (sweeps < max_sweeps && !active.empty()) {
      double amax = 0.0;
      for (const int j : active) {
        amax = std::max(amax, update_coord(j));
      }
      ++sweeps;
      if (amax < inner_tol) break;
    }
    // Loop back to a full sweep to check for support violations.
  }
  if (!converged) kkt = kkt_residual();

  return Rcpp::List::create(Rcpp::Named("x") = x,
                            Rcpp::Named("sweeps") = sweeps,
                            Rcpp::Named("kkt") = kkt,
                            Rcpp::Named("converged") = converged);
}
