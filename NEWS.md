# genpca 0.2.0

### Release review fixes

* Randomized sketches use Gaussian draws and use a deterministic complete basis
  when the sketch spans a full matrix dimension. Metric orthonormalization
  verifies the original Gram matrix, corrects residual normalization error,
  and removes dependent sketch directions without inventing rank through
  jitter. This fixes seed- and C++-library-dependent component loss on Windows.
* `gpca_mle()` now separates the penalty effect of reciprocal rescaling
  (`loglik_rescale_delta`, exactly zero with `scale_fix = "none"`) from final
  refitting and objective reevaluation (`loglik_refit_delta`). The reported
  likelihood remains evaluated at the returned fit and metrics.

* Strictly positive diagonal weights are retained in both forward and inverse
  metric factors. Previously a small positive weight could contribute a large
  singular value while its recovered loading was zero. The fix covers GPCA,
  covariance PCA and PLS operators, including base-matrix diagonal inputs.
* Both deflation implementations now separate `rank_rtol` from the iteration
  `threshold`. All backends also apply a common final component filter. Empty
  C++ deflation results return empty singular-value and variance vectors.
* Incomplete eigencore results raise `genpca_solver_nonconvergence`. Existing
  dense fallbacks handle this error; operator callers without a fallback stop
  instead of returning an unchecked fit.
* Explicit `"clip"` repairs remove negative eigenvalues even within the PSD
  validation tolerance and report the change, up to reconstruction roundoff.
* `geigen_cov()` documents its projected equation for singular metrics and the
  range constraint on its vectors. Commuting metrics and covariances need not
  give the same leading component under GMD and generalized eigenanalysis.
* Rank-cutoff equations render in both PDF and HTML manuals. Local diagnostic
  artifacts are excluded from source packages.

### Numerical semantics (audit of 2026-09)

* All rank and validity decisions are now relative to the scale of the
  problem, so results are invariant to rescaling `X` or a metric. New
  `rank_rtol` argument in `genpca()`, `genpca_cov()` and the internal
  solvers: component `j` is dropped when `d_j <= rank_rtol * d_1` (default
  `1e-6`), for every method. Metric validation and null-space detection use
  a separate relative tolerance, `sqrt(.Machine$double.eps)`. The previous
  absolute cutoffs (`1e-8` on eigenvalues, `1e-9` on singular values, a
  `1e-8` total-variance floor, an absolute `jitter_metric`) could change the
  number of components returned, or return wrong values from the randomized
  backend, when `X` was rescaled.
* Positive semi-definiteness and positive definiteness are now distinct
  checks (`is_psd()`, `is_pd()`, both shifted Cholesky probes with relative
  tolerance). `ensure_spd()` always returns a positive definite matrix and
  forwards its tolerance; previously a zero or singular matrix could be
  returned unchanged.
* Asymmetric metrics are an error under every `constraints_remedy`. The
  relative asymmetry `||A - A'|| / ||A||` is measured; below `1e-10` the two
  triangles are averaged, above it the call stops. Previously one triangle
  was silently used, so the result depended on which triangle the caller
  had filled. The same applies to `C` and `R` in `genpca_cov()`.
* `constraints_remedy = "clip"` now guarantees a PSD result; the previous
  tolerant fast path could return a matrix with a small negative eigenvalue.
* Diagonal weights (vectors or diagonal matrices) are validated identically
  by every backend: entries within the relative tolerance of zero are set to
  exactly zero with a message; more negative entries are subject to
  `constraints_remedy`. Previously the eigen path rejected what the other
  backends silently accepted or clamped.
* `method = "eigen"` never truncates a metric any more. A positive definite
  general metric is factored exactly by Cholesky at any size; a singular
  general metric needs a dense eigendecomposition, which is refused above
  `maxeig` rows (default now `5000`) with a message naming the
  alternatives. The previous behaviour, keeping only the `maxeig` largest
  eigenvectors of the metric before looking at `X`, could discard the
  dominant component entirely (a direction with a small metric eigenvalue
  and large data variance). `warn_approx` is deprecated and ignored.
* `genpca_cov(method = "gmd")` requires `C` to be PSD within tolerance and
  no longer clamps negative metric eigenvalues silently. Its `tol` argument
  is deprecated in favor of `rank_rtol` and `metric_rtol`.
* New exported `repair_metric()` returns a repaired metric together with a
  `metric_repair_report` (minimum eigenvalue before and after, Gershgorin
  bound, applied shift, rank, condition number).

### Explicit metric repair (breaking)

* `genpca()` now defaults to `constraints_remedy = "error"`: an indefinite
  metric stops the fit instead of being silently shifted. Opt into a repair
  with `constraints_remedy = "ridge"`, `"clip"` or `"identity"`; every
  repair that changes the metric emits a warning of class
  `genpca_metric_repaired` whose `report` field is the `repair_metric()`
  diagnostic. `gpca_mle()` inherits the new default.
* `genpls()`, `genplsc()` and `gplssvd_op()` gain the same
  `constraints_remedy` argument (default `"error"`). Previously their
  metric operators applied a Gershgorin ridge to any indefinite weight
  matrix with no message and no way to opt out.
* `verbose = TRUE` is no longer needed to find out that a metric was
  replaced by the identity; the warning always fires.

### Covariance interface and `gpca_mle()`

* `genpca_cov()` is now the GMD estimator only (eigendecomposition of
  `R^{1/2} C R^{1/2}`). The generalized eigenproblem `C v = lambda R v`, a
  different estimator, is the new exported
  `geigen_cov()`. `genpca_cov(method = "geigen")` still works for one
  release and forwards with a deprecation warning.
* `gpca_mle()` defaults to `scale_fix = "none"`. The ridge penalty already
  identifies the overall scale of the learned metrics, and the previous
  default `"trace"` rescaled the metrics after convergence in a way the
  penalized objective is not invariant to, while still reporting the
  pre-rescale value. `loglik` is now always the penalized objective at the
  returned metrics; `loglik_unpenalized` and `loglik_rescale_delta` (the
  change caused by a `"trace"`/`"det"` rescale) are new.

### Solver backend

* The iterative eigen/SVD backend is now the **eigencore** package;
  **RSpectra** is no longer a dependency. `method = "spectra"` computes the
  top-k singular triplets of the metric-whitened data `F_M' X F_A`
  (`M = F_M F_M'`, `A = F_A F_A'`) as an implicit operator: diagonal, dense
  Cholesky, sparse (CHOLMOD) Cholesky and eigen factors (for singular
  metrics) are all supported, and a dense SVD is used when the iterative
  solver does not converge. The C++ Spectra kernel and the `LinkingTo:
  RSpectra` dependency are gone. `gmd_fast_cpp()` remains as an alias of the
  new `gmd_spectra()`; its `maxit` and `seed` arguments are ignored.
* `genpls()`, `genplsc()` and `gplssvd_op()` accept
  `svd_backend = "eigencore"` (new default); `"RSpectra"` is accepted as a
  deprecated alias.
* Iterative solves never disturb the caller's `.Random.seed`.

# genpca 0.1.0

* Initial CRAN submission.
* Generalized PCA with row metric `M` and column metric `A`
  (`genpca()`), following Allen, Grosenick & Taylor (2014).
* Covariance-based GPCA from a precomputed `C = X' M X`
  (`genpca_cov()`), with `eigen` and `gmd` paths.
* Multiple computational backends: `eigen`, `spectra` (matrix-free C++
  via RSpectra), `randomized`, and `deflation`. The `auto` heuristic
  picks among them.
* Generalized PLS / PLS-SVD on two blocks (`genpls()`, `genplsc()`)
  and an operator-level interface (`gplssvd_op()`) that avoids
  materialising whitened matrices.
* Sparse functional PCA (`sfpca()`), regularised PLS (`rpls()`), and
  matrix-normal PCA via maximum residual likelihood (`mnpca_mrl()`).
  `sfpca()` solves each rank-1 subproblem exactly via C++ coordinate
  descent in the constraint form of Allen & Weylandt (2019), deflates
  implicitly (sparse inputs never densify), and selects sparsity
  penalties per component by BIC along a warm-started `lambda` path;
  smoothness weights default to the scale-free `1 / lambda_max(Omega)`.
  The former `uthresh`/`vthresh` quantile heuristics are deprecated.
  `sfpca()` now returns a `multivarious` `bi_projector` (class `"sfpca"`),
  so `scores()`, `components()`, `sdev()`, and `reconstruct()` work as
  they do for `genpca()`. The pre-0.1 list fields `$d` and `$u` remain
  readable but emit a deprecation warning.
* Maximum-likelihood metric learning (`gpca_mle()`).
* SPD constraint handling with `"ridge"`, `"clip"`, and `"identity"`
  remediation strategies.
* Vignettes covering getting started, metric recipes, scaling, and a
  reference implementation for GPLSSVD.

### Bug fixes and behavior changes

* `gpca_mle()` now optimizes a single penalized (MAP) matrix-normal
  objective with exact block updates (sequential flip-flop covariance
  updates; the ridge `lambda` is part of the objective), so
  `loglik_path` is monotone non-decreasing. `scale_fix` is applied once
  at exit as a joint, likelihood-preserving rescale (row covariance
  normalized, factor absorbed into the column covariance) instead of
  normalizing both factors independently during iteration, and the
  returned `fit` is computed with the returned canonicalized metrics.
  Reported log-likelihood values now include the penalty term, so they
  differ numerically from previous versions.
* `scores()` on `genpca` fits now matches Allen, Grosenick & Taylor (2014)'s
  convention (`z_k = X A ov_k = ou_k d_k`) and is identical to
  `project(fit, X)` on the training data; previously the two could disagree.
* `reconstruct.genpca(fit, colind = )` now applies the inverse
  pre-processing transform to the selected columns, rather than to the
  first `length(colind)` columns.
* Deflation convergence and degeneracy thresholds (the `threshold` argument
  of `genpca(method = "deflation")` and the internal SFPCA solver) are now
  relative to the scale of the problem rather than absolute, so results are
  invariant to rescaling the input data.
* `constraints_remedy = "clip"` now performs a real spectral clip to the
  PSD cone (previously it did not clip); it requires a dense
  eigendecomposition and refuses sparse input larger than 2000 rows/columns,
  recommending `"ridge"` instead.
* Fixed PSD validation for dense (non-diagonal) constraint matrices so that
  indefinite input is reliably rejected/repaired rather than silently
  accepted.
* `genpca(method = "randomized")` no longer mutates the caller's
  `.Random.seed`; the C++ kernel seeds its own stream from
  `seed_randomized`, which now fully controls reproducibility.
* Added input validation (dimensions, finiteness, penalty/`scad_a` range,
  strictly positive diagonal) to the internal C++ coordinate-descent solver
  used by `sfpca()`.
* Fixed a cache-key collision in the internal matrix-decomposition cache:
  the digest previously used rounded matrix entries, which could collide
  for distinct matrices (e.g. a metric and its ridge-remediated version)
  and return the wrong cached factor; the digest now hashes exact bytes.
* `genpls()`'s `$ncomp` field now reports the number of components actually
  extracted, which may be less than the requested `ncomp`.
