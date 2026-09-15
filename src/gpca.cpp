#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Generalized power-iteration deflation (Allen, Grosenick & Taylor, 2014):
//   uhat = Xhat R v;  u = uhat / sqrt(u' Q u)
//   vhat = Xhat' Q u; v = vhat / sqrt(v' R v)
//   d_i  = u' Q X R v
// Degenerate directions (zero norm in the metric) or near-zero singular
// values stop extraction early and results are trimmed, matching the R
// implementation gmd_deflationR().
//
// Stopping criteria are RELATIVE to the scale of X in the (Q, R) metric
// (scale_ref = ||X||_{Q,R}), so results are invariant to rescaling X.
//
// Warnings are collected into the returned list ("warnings") and emitted by
// the R wrapper: calling Rf_warning() from C++ longjmps past Armadillo
// destructors under options(warn = 2), leaking every live allocation.
//
// The residual is applied implicitly.  After j components have been extracted:
//   X_j  w = X w  - U_j (d_j * (V_j' w))
//   X_j' z = X' z - V_j (d_j * (U_j' z))
// This keeps sparse X sparse and avoids a full dense residual copy.

// large dimension should be in ROWS (if columns is X, then pass (t(X), R,))

template <typename MatX>
List gmd_deflation_impl(const MatX &X, const arma::sp_mat &Q, const arma::sp_mat &R,
                        int k, double thr=1e-7, int maxit=500, bool verbose=false, double rank_rtol=1e-6) {
  if (maxit < 1) {
    stop("maxit must be >= 1.");
  }
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat ugmd(n, k, fill::zeros);
  arma::mat vgmd(p, k, fill::zeros);
  arma::vec dgmd(k, fill::zeros);
  arma::vec propv(k, fill::zeros);
  arma::vec cumv(k, fill::zeros);
  std::vector<std::string> warns;

  // qrnorm = tr(X' Q X R) = accu((Q X) % (X R)); avoids any n x n / p x p
  // dense intermediate (trace is cyclic, so one expression covers both shapes).
  double qrnorm = arma::accu((Q * X) % (X * R));
  double scale_ref = 1.0;
  if (!std::isfinite(qrnorm) || qrnorm <= 0.0) {
    warns.push_back("Total generalized variance is not positive; explained variance may be unstable.");
    qrnorm = 1.0;
  } else {
    scale_ref = std::sqrt(qrnorm);
  }

  int k_found = 0;

  auto residual_mv = [&](const arma::vec& w, const int n_prev) -> arma::vec {
    arma::vec y = X * w;
    if (n_prev > 0) {
      const arma::uword nprev = static_cast<arma::uword>(n_prev);
      const arma::uword last = nprev - 1;
      arma::vec coeff = dgmd.head(nprev) % (vgmd.cols(0, last).t() * w);
      y -= ugmd.cols(0, last) * coeff;
    }
    return y;
  };

  auto residual_t_mv = [&](const arma::vec& z, const int n_prev) -> arma::vec {
    arma::vec y = X.t() * z;
    if (n_prev > 0) {
      const arma::uword nprev = static_cast<arma::uword>(n_prev);
      const arma::uword last = nprev - 1;
      arma::vec coeff = dgmd.head(nprev) % (ugmd.cols(0, last).t() * z);
      y -= vgmd.cols(0, last) * coeff;
    }
    return y;
  };

  const double norm_floor = std::numeric_limits<double>::epsilon() * scale_ref;

  for (int i=0; i<k; i++) {
    Rcpp::checkUserInterrupt();

    // Fresh random start for every component: the previous component's
    // converged vectors are (by construction) annihilated by the deflated
    // operator, so reusing them would amplify roundoff noise into the
    // starting direction.
    arma::vec u = randn(n);
    arma::vec v = randn(p);

    double err = 1;
    int iter = 0;
    bool degenerate = false;
    while (err > thr && iter < maxit) {
      iter += 1;
      arma::vec oldu = u;
      arma::vec oldv = v;

      arma::vec uhat = residual_mv(R * v, k_found);
      double u_norm = std::sqrt(arma::as_scalar(uhat.t() * (Q * uhat)));
      if (!std::isfinite(u_norm) || u_norm <= norm_floor) { degenerate = true; break; }
      u = uhat / u_norm;

      arma::vec vhat = residual_t_mv(Q * u, k_found);
      double v_norm = std::sqrt(arma::as_scalar(vhat.t() * (R * vhat)));
      if (!std::isfinite(v_norm) || v_norm <= norm_floor) { degenerate = true; break; }
      v = vhat / v_norm;

      err = arma::accu(arma::square(oldu - u)) + arma::accu(arma::square(oldv - v));
    }

    if (degenerate) {
      if (verbose) {
        Rcout << "Degenerate direction at component " << (i + 1) << "; stopping deflation." << std::endl;
      }
      warns.push_back("Deflation stopped early at component " + std::to_string(i + 1) +
                      " (degenerate direction).");
      break;
    }

    if (iter >= maxit && err > thr) {
      if (verbose) {
        Rcout << "Power iteration reached maxit for component " << (i + 1) << "." << std::endl;
      }
      warns.push_back("Power iteration reached maxit at component " + std::to_string(i + 1) +
                      " (maxit=" + std::to_string(maxit) + ").");
    }

    double d_i = arma::as_scalar((Q * u).t() * residual_mv(R * v, k_found));
    // Relative cutoff: stop once the residual singular value is negligible
    // compared to the largest one extracted; the first establishes the scale.
    double d_ref = (k_found > 0) ? std::fabs(dgmd(0)) : std::fabs(d_i);
    if (!std::isfinite(d_i) || d_i <= 0.0 || d_i <= rank_rtol * d_ref) {
      warns.push_back("Deflation stopped early at component " + std::to_string(i + 1) +
                      " (singular value near zero).");
      break;
    }

    dgmd(k_found) = d_i;
    ugmd.col(k_found) = u;
    vgmd.col(k_found) = v;
    propv(k_found) = d_i * d_i / qrnorm;
    cumv(k_found) = (k_found > 0) ? (cumv(k_found - 1) + propv(k_found)) : propv(k_found);
    k_found += 1;
  }

  return List::create(
    Named("d") = arma::vec(dgmd.head(k_found)),
    Named("v") = vgmd.head_cols(k_found),
    Named("u") = ugmd.head_cols(k_found),
    Named("k") = k_found,
    Named("cumv") = arma::vec(cumv.head(k_found)),
    Named("propv") = arma::vec(propv.head(k_found)),
    Named("warnings") = wrap(warns)
  );
}

//[[Rcpp::export]]
List gmd_deflation_cpp(const arma::mat &X, const arma::sp_mat &Q, const arma::sp_mat &R,
                       int k, double thr=1e-7, int maxit=500, bool verbose=false, double rank_rtol=1e-6) {
  return gmd_deflation_impl(X, Q, R, k, thr, maxit, verbose, rank_rtol);
}

//[[Rcpp::export]]
List gmd_deflation_cpp_sp(const arma::sp_mat &X, const arma::sp_mat &Q, const arma::sp_mat &R,
                          int k, double thr=1e-7, int maxit=500, bool verbose=false, double rank_rtol=1e-6) {
  return gmd_deflation_impl(X, Q, R, k, thr, maxit, verbose, rank_rtol);
}
