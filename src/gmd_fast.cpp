// Copyright (c) 2025 genpca contributors
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <random>
// [[Rcpp::depends(RcppArmadillo)]]

// ---- helpers ---------------------------------------------------------------

static inline arma::mat random_normal_matrix(const arma::uword nrow,
                                           const arma::uword ncol,
                                           const unsigned int seed) {
  std::mt19937 gen(seed);
  std::normal_distribution<double> normal(0.0, 1.0);
  arma::mat out(nrow, ncol);
  for (arma::uword j = 0; j < ncol; ++j) {
    for (arma::uword i = 0; i < nrow; ++i) {
      out(i, j) = normal(gen);
    }
  }
  return out;
}

template <typename MatM>
static arma::mat metric_orthonormalize_cpp(const arma::mat& A,
                                           const MatM& M,
                                           const double jitter,
                                           const double tol) {
  if (A.n_cols == 0) {
    return arma::mat(A.n_rows, 0, arma::fill::zeros);
  }

  arma::mat MA = M * A;
  arma::mat G = A.t() * MA;
  G = 0.5 * (G + G.t());

  // A jittered Cholesky is a candidate preconditioner, not proof that the
  // original metric Gram is full rank. Certify against M before accepting it.
  const double eps = std::numeric_limits<double>::epsilon();
  const double gscale = G.n_rows ? std::max(0.0, G.diag().max()) : 0.0;
  arma::mat Greg = G;
  Greg.diag() += std::max(0.0, jitter) * gscale;
  arma::mat C;
  if (arma::chol(C, Greg, "upper")) {
    arma::mat Xt = arma::solve(arma::trimatl(C.t()), A.t(), arma::solve_opts::fast);
    arma::mat B = Xt.t();
    arma::mat check = B.t() * (M * B);
    if (check.is_finite() &&
        arma::norm(check - arma::eye(check.n_rows, check.n_cols), "inf") <= 1e-8) {
      // Remove the small normalization error left by jitter (CholeskyQR2).
      arma::mat correction;
      if (arma::chol(correction, check, "upper")) {
        arma::mat corrected = arma::solve(arma::trimatl(correction.t()), B.t(),
                                          arma::solve_opts::fast);
        return corrected.t();
      }
    }
  }

  // Remove linear dependence in Euclidean coordinates before forming a
  // metric Gram. Column scaling and SVD avoid squaring the raw sketch's
  // condition number. The user singular-value cutoff belongs to the final
  // decomposition, not to this internal basis or the metric spectrum.
  arma::mat scaled = A;
  for (arma::uword j = 0; j < scaled.n_cols; ++j) {
    const double scale = arma::norm(scaled.col(j), 2);
    if (scale > 0.0) scaled.col(j) /= scale;
  }
  arma::mat basis, right;
  arma::vec singular;
  if (!arma::svd_econ(basis, singular, right, scaled))
    Rcpp::stop("SVD failed while orthonormalizing the randomized sketch.");
  const double floor = eps * std::max(A.n_rows, A.n_cols) *
    (singular.n_elem ? singular.max() : 0.0);
  arma::uvec independent = arma::find(singular > floor && singular > 0.0);
  if (independent.n_elem == 0) return arma::mat(A.n_rows, 0, arma::fill::zeros);
  basis = basis.cols(independent);
  G = basis.t() * (M * basis);
  G = 0.5 * (G + G.t());
  arma::vec eval;
  arma::mat evec;
  if (!arma::eig_sym(eval, evec, G))
    Rcpp::stop("Metric eigendecomposition failed while orthonormalizing the sketch.");
  const double emax = eval.n_elem ? std::max(0.0, eval.max()) : 0.0;
  arma::uvec keep = arma::find(eval > eps * G.n_rows * emax && eval > 0.0);
  if (keep.n_elem == 0) return arma::mat(A.n_rows, 0, arma::fill::zeros);
  return basis * evec.cols(keep) * arma::diagmat(1.0 / arma::sqrt(eval.elem(keep)));
}

template <typename MatQ, typename MatR>
static void randomized_polish_cpp(const arma::mat& X,
                                  const MatQ& Q,
                                  const MatR& R,
                                  arma::mat& U,
                                  arma::mat& V,
                                  arma::vec& d,
                                  const int n_polish,
                                  const double jitter,
                                  const double tol,
                                  const double polish_tol) {
  if (n_polish <= 0 || U.n_cols == 0) return;

  arma::vec d_prev;
  for (int it = 0; it < n_polish; ++it) {
    arma::mat Y = X * (R * V);
    U = metric_orthonormalize_cpp(Y, Q, jitter, tol);
    arma::mat Z = X.t() * (Q * U);
    V = metric_orthonormalize_cpp(Z, R, jitter, tol);
    arma::mat RV = R * V;
    // T = U^T Q X R V; reuse Z = X^T Q U to avoid another pass through X.
    arma::mat T = Z.t() * RV;

    arma::mat P, Vrot;
    arma::vec s;
    if (!arma::svd(P, s, Vrot, T)) {
      Rcpp::stop("SVD failed in randomized polish.");
    }
    U = U * P;
    V = V * Vrot;
    d = s;

    if (polish_tol > 0.0 && d_prev.n_elem == d.n_elem && d.n_elem > 0) {
      arma::vec denom = arma::max(arma::abs(d_prev), arma::vec(d_prev.n_elem, arma::fill::value(1e-12)));
      arma::vec rel = arma::abs(d - d_prev) / denom;
      double max_rel = rel.max();
      if (std::isfinite(max_rel) && max_rel < polish_tol) {
        break;
      }
    }
    d_prev = d;
  }
}

template <typename MatQ, typename MatR>
Rcpp::List gmd_randomized_impl(const arma::mat& X,
                               const MatQ& Q,
                               const MatR& R,
                               const int k,
                               const int oversample,
                               const int n_power,
                               const int n_polish,
                               const double jitter,
                               const double tol,
                               const double polish_tol,
                               const unsigned int seed) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  const int k_use = std::max(0, std::min(k, std::min(n, p)));

  if (k_use == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("k") = 0
    );
  }

  const int ell = std::max(1, std::min(std::min(n, p), k_use + std::max(0, oversample)));
  arma::mat U0;
  if (ell == n) {
    // The sketch covers the entire row space: use it directly, without a
    // random draw that can accidentally lose a direction.
    U0 = metric_orthonormalize_cpp(arma::eye(n, n), Q, jitter, tol);
  } else {
    arma::mat Omega = ell == p ? arma::eye(p, p) :
      random_normal_matrix(static_cast<arma::uword>(p),
                           static_cast<arma::uword>(ell), seed);
    arma::mat Y = X * (R * Omega);
    for (int it = 0; it < n_power; ++it) {
      arma::mat Utmp = metric_orthonormalize_cpp(Y, Q, jitter, tol);
      if (Utmp.n_cols == 0) break;
      arma::mat Z = X.t() * (Q * Utmp);
      Y = X * (R * Z);
    }
    U0 = metric_orthonormalize_cpp(Y, Q, jitter, tol);
  }
  if (U0.n_cols == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("k") = 0
    );
  }

  arma::mat B = X.t() * (Q * U0);
  arma::mat RB = R * B;
  arma::mat G = B.t() * RB;
  G = 0.5 * (G + G.t());

  arma::vec eval;
  arma::mat Sfull;
  if (!arma::eig_sym(eval, Sfull, G)) {
    Rcpp::stop("eig_sym failed in randomized solver.");
  }
  arma::uvec ord = arma::sort_index(eval, "descend");
  const arma::uword kk = std::min<arma::uword>(static_cast<arma::uword>(k_use), ord.n_elem);
  if (kk == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("k") = 0
    );
  }
  arma::uvec take = ord.head(kk);

  arma::vec lam = eval.elem(take);
  for (arma::uword i = 0; i < lam.n_elem; ++i) {
    if (lam(i) < 0.0) lam(i) = 0.0;
  }
  arma::vec d = arma::sqrt(lam);
  arma::mat S = Sfull.cols(take);

  arma::mat U = U0 * S;
  arma::mat V = B * S;
  const double dmax0 = d.n_elem ? d.max() : 0.0;
  for (arma::uword i = 0; i < d.n_elem; ++i) {
    if (d(i) > 0.0 && d(i) > tol * dmax0) {
      V.col(i) /= d(i);
    } else {
      V.col(i).zeros();
    }
  }

  randomized_polish_cpp(X, Q, R, U, V, d, n_polish, jitter, tol, polish_tol);

  const double dmax = d.n_elem ? d.max() : 0.0;
  arma::uvec keep = arma::find(d > 0.0 && d > tol * dmax);
  if (keep.n_elem == 0) {
    return Rcpp::List::create(
      Rcpp::Named("u") = arma::mat(X.n_rows, 0, arma::fill::zeros),
      Rcpp::Named("v") = arma::mat(X.n_cols, 0, arma::fill::zeros),
      Rcpp::Named("d") = arma::vec(),
      Rcpp::Named("k") = 0
    );
  }

  U = U.cols(keep);
  V = V.cols(keep);
  d = d.elem(keep);

  return Rcpp::List::create(
    Rcpp::Named("u") = U,
    Rcpp::Named("v") = V,
    Rcpp::Named("d") = d,
    Rcpp::Named("k") = static_cast<int>(d.n_elem)
  );
}

// ---- Exported entry points (randomized backend) ----------------------------

// [[Rcpp::export]]
Rcpp::List gmd_randomized_cpp_dn(const arma::mat& X,
                                 const arma::mat& Q,
                                 const arma::mat& R,
                                 const int k,
                                 const int oversample = 20,
                                 const int n_power = 1,
                                 const int n_polish = 0,
                                 const double jitter = 1e-10,
                                 const double tol = 1e-9,
                                 const double polish_tol = 0.0,
                                 const int seed = 1234) {
  return gmd_randomized_impl(
    X, Q, R, k, oversample, n_power, n_polish, jitter, tol, polish_tol, static_cast<unsigned int>(seed)
  );
}

// [[Rcpp::export]]
Rcpp::List gmd_randomized_cpp_sp(const arma::mat& X,
                                 const arma::sp_mat& Q,
                                 const arma::sp_mat& R,
                                 const int k,
                                 const int oversample = 20,
                                 const int n_power = 1,
                                 const int n_polish = 0,
                                 const double jitter = 1e-10,
                                 const double tol = 1e-9,
                                 const double polish_tol = 0.0,
                                 const int seed = 1234) {
  return gmd_randomized_impl(
    X, Q, R, k, oversample, n_power, n_polish, jitter, tol, polish_tol, static_cast<unsigned int>(seed)
  );
}

// [[Rcpp::export]]
Rcpp::List gmd_randomized_cpp_qsp_rdn(const arma::mat& X,
                                      const arma::sp_mat& Q,
                                      const arma::mat& R,
                                      const int k,
                                      const int oversample = 20,
                                      const int n_power = 1,
                                      const int n_polish = 0,
                                      const double jitter = 1e-10,
                                      const double tol = 1e-9,
                                      const double polish_tol = 0.0,
                                      const int seed = 1234) {
  return gmd_randomized_impl(
    X, Q, R, k, oversample, n_power, n_polish, jitter, tol, polish_tol, static_cast<unsigned int>(seed)
  );
}

// [[Rcpp::export]]
Rcpp::List gmd_randomized_cpp_qdn_rsp(const arma::mat& X,
                                      const arma::mat& Q,
                                      const arma::sp_mat& R,
                                      const int k,
                                      const int oversample = 20,
                                      const int n_power = 1,
                                      const int n_polish = 0,
                                      const double jitter = 1e-10,
                                      const double tol = 1e-9,
                                      const double polish_tol = 0.0,
                                      const int seed = 1234) {
  return gmd_randomized_impl(
    X, Q, R, k, oversample, n_power, n_polish, jitter, tol, polish_tol, static_cast<unsigned int>(seed)
  );
}
