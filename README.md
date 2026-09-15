# genpca

An R package for **Generalized Principal Component Analysis (GPCA)** and related matrix decompositions when data is observed in non-Euclidean inner-product spaces. 

## What is Generalized PCA?

Standard PCA assumes all observations and variables are equally important and that Euclidean distance is the appropriate similarity measure. However, many real-world datasets violate these assumptions:

- **Weighted observations**: Survey data where rows represent different population sizes
- **Variable precision**: Measurements with different accuracies or importance
- **Correlated features**: Spatial/temporal data with known dependency structures  
- **Domain-specific geometry**: Functional data, shape analysis, or other specialized metrics

GPCA extends standard PCA by incorporating row and column metrics (M and A) that encode prior knowledge about data structure, following the framework of Allen, Grosenick & Taylor (2014). The package also implements generalized PLS methods building on ideas from Beaton et al. (2016).

## Key Features

### Core Functionality
- **Generalized PCA** (`genpca`): Decomposition with row metric M and column metric A
- **Covariance-based GPCA** (`genpca_cov`): Direct analysis of pre-computed covariance matrices
- **Multiple computational backends**:
  - `eigen`: Direct eigendecomposition for small-to-medium problems
  - `spectra`: Iterative partial SVD of the metric-whitened data (via eigencore) for large problems
  - `randomized`: Block-sketch approximation for wide low-rank settings
  - `deflation`: Sequential extraction for memory-constrained scenarios

### Advanced Methods
- **Generalized PLS/PLS-SVD** (`genpls`/`genplsc`): Two-block analysis with metrics
- **Operator-level computations** (`gplssvd_op`): Efficient PLS without materializing whitened matrices
- **Constraint handling**: Strict, scale-relative validation of metric matrices (symmetry and positive semi-definiteness); repairs are explicit and reported (`repair_metric()`, `constraints_remedy`)

### Related Decompositions
- **Sparse functional PCA** (`sfpca`): Rank-1 components with sparsity and spatial smoothness penalties, e.g. `sfpca(X, K = 2, spat_cds = coords)`. See `?sfpca`.
- **Regularized PLS** (`rpls`): Single-metric penalized PLS (Allen et al. 2013) with `"l1"` or `"ridge"` penalties, e.g. `rpls(X, Y, K = 3, lambda = 0.1, penalty = "l1")`. See `?rpls`.
- **Matrix-normal PCA** (`mnpca_mrl`): Low-rank factorization with sparse row/column precision matrices under matrix-normal noise, e.g. `mnpca_mrl(Y, ncomp = 3, lambda_row = 0.05, lambda_col = 0.05)`. See `?mnpca_mrl`.
- **Metric learning** (`gpca_mle`): Experimental alternation of GPCA with maximum-likelihood estimation of the row/column metrics themselves, e.g. `gpca_mle(X, ncomp = 2)`. See `?gpca_mle`.

### Integration
- Full compatibility with the `multivarious` package ecosystem
- Unified interface for preprocessing, projection, reconstruction, and transfer learning
- Support for sparse matrices via the `Matrix` package

## Backend Selection Guide

The default backend is `method = "eigen"`. `"auto"` is opt-in: pass it
explicitly to let a heuristic pick among `"eigen"`, `"spectra"`, and
`"randomized"` based on problem shape and constraint structure.

| Method | Best for | Pros | Cons |
|---|---|---|---|
| `eigen` | Small/medium problems, exact reference runs | Most stable reference behavior | Builds larger intermediate matrices; for very large sparse constraints may rely on truncated eigensolve (`maxeig`) and become slower/approximate |
| `spectra` | Large problems with few components | Iterative partial SVD of the whitened operator (eigencore); lower memory | Both metrics are factored once (Cholesky or sparse Cholesky); a large dense general metric costs one factorization |
| `randomized` | Wide (`p >> n`), sparse-metric, low-rank workloads | Often fastest in wide settings; block GEMM/SpMM path | Approximate by design; tune `oversample`, `n_power`, `n_polish` |
| `deflation` | Few components with limited memory | Low memory, component-by-component extraction | Can converge slowly; C++ path currently expects sparse metrics |
| `auto` | Default production usage | Chooses among `eigen`/`spectra`/`randomized` heuristically | Heuristics may not be optimal for every hardware/data regime |

## Installation

```r
# install.packages("devtools")
devtools::install_github("bbuchsbaum/genpca")
```

You’ll also want these runtime dependencies installed:

```r
install.packages(c("Matrix", "eigencore", "multivarious"))
# Optional for some utilities / tests
install.packages(c("irlba", "knitr", "rmarkdown"))
```

## Quick Start

### Basic Usage

```r
library(genpca)
set.seed(1)
X <- matrix(rnorm(200 * 50), 200, 50)

# Standard PCA (identity metrics by default)
fit <- genpca(X, ncomp = 5, preproc = multivarious::center())
fit$sdev                             # generalized singular values;
                                      # with identity metrics this is
                                      # prcomp(X)$sdev * sqrt(nrow(X) - 1)
head(multivarious::scores(fit))      # scores (n × k)
head(multivarious::components(fit))  # loadings (p × k)
```

### Weighted GPCA

```r
# Example: Survey data with population weights
library(Matrix)
pop_weights <- runif(nrow(X), 0.5, 1.5)  # population sizes
M <- Diagonal(nrow(X), x = pop_weights)  # row metric

# Variable importance weights  
var_importance <- c(rep(2, 10), rep(1, 30), rep(0.5, 10))
A <- Diagonal(ncol(X), x = var_importance)  # column metric

fit_weighted <- genpca(X, M = M, A = A, ncomp = 5, 
                      preproc = multivarious::center())
```

### Covariance-based GPCA

```r
# When you have pre-computed covariance C = X'MX. fit_weighted above used
# preproc = multivarious::center(), so C must be built from centered X for
# the two fits to agree exactly.
Xc <- scale(X, center = TRUE, scale = FALSE)
C  <- crossprod(Xc, M %*% Xc)
fit_cov <- genpca_cov(C, R = A, ncomp = 5, method = "gmd")
# Mathematically equivalent to fit_weighted above
all.equal(fit_weighted$sdev, fit_cov$d, tolerance = 1e-8)
```

### Generalized PLS

```r
# Two-block analysis with canonical PLS
Y <- matrix(rnorm(200 * 20), 200, 20)
pls <- genpls(X, Y, ncomp = 3, 
              preproc_x = multivarious::center(), 
              preproc_y = multivarious::center())
pls$d                    # singular values of the whitened cross-product
                          # Xe'Ye (covariance scale, NOT canonical correlations)
dim(pls$vx); dim(pls$vy) # X/Y weight matrices

# To get canonical-correlation-like quantities, correlate the latent
# variable pairs directly:
diag(cor(pls$lx, pls$ly))
```

## Understanding Metrics in GPCA

The metrics M and A define inner products and distances in the observation and variable spaces:

### Row Metric M (n × n)
- Defines relationships between observations
- Inner product: `||x||_M^2 = x^T M x`
- Distance: `d_M(x,y)^2 = (x−y)^T M (x−y)`
- **Common choices:**
  - Identity: Standard equal weighting
  - Diagonal: Population/sample weights
  - Precision matrix: Account for observation correlations
  - Kernel matrices: Encode similarity structures

### Column Metric A (p × p)  
- Defines relationships between variables
- Inner product: `||v||_A^2 = v^T A v`
- Distance: `d_A(v,w)^2 = (v−w)^T A (v−w)`
- **Common choices:**
  - Identity: Standard equal importance
  - Diagonal: Variable weights/importance
  - Covariance/precision: Variable dependencies
  - Graph Laplacian: Spatial/temporal smoothness

When `M = I` and `A = I`, GPCA reduces to standard PCA.

## Documentation

- Vignettes
  - Getting Started with genpca
  - GPCA Metrics: Building M and A
  - Modelling Structured Noise
  - GPCA at Scale and Special Cases
  - Generalized PLS-SVD: Explicit Whitening Reference

Build locally:

```r
devtools::build_vignettes()
browseVignettes("genpca")
```

## Testing and guarantees

- Eigen vs spectra (eigencore): unit tests assert tight agreement on modest problems (sdev within 1e‑6, scores within 1e‑5 up to sign).
- Deflation vs Eigen: additional tests (n≈60, p≈40, k=8) assert:
  - sdev within 1e‑4,
  - subspace agreement via principal angles,
  - scores/components close after Procrustes alignment (≈ 1e‑3 relative).

Run tests locally:

```r
library(testthat)
library(pkgload)
pkgload::load_all()
testthat::test_dir("tests/testthat")
```

## References

The methods in this package are based on:

- **Allen, G. I., Grosenick, L., & Taylor, J. (2014).** A generalized least-square matrix decomposition. *Journal of the American Statistical Association*, 109(505), 145-159. doi:10.1080/01621459.2013.852978

- **Beaton, D., ADNI, et al. (2016).** Generalized partial least squares: A framework for simultaneously capturing common and individual variation. *NeuroImage*, 141, 346-363. doi:10.1016/j.neuroimage.2016.07.034

For additional theoretical background on generalized decompositions, see:

- **Beaton, D. (2020).** Generalized eigen, singular value, and partial least squares decompositions: The GSVD package. *arXiv preprint* arXiv:2010.14734.

## License

MIT (see LICENSE).

## Contributing

Issues and PRs welcome. Please open a ticket with a minimal example, your R session info, and (if relevant) a pointer to the metric matrices that reproduce the behavior.



<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The defaults are configured via `params$family` and `params$preset` (family = 'red', preset = 'interaction'). The pkgdown site uses `template: { package: albersdown }` together with generated `pkgdown/extra.css` and `pkgdown/extra.js` so the theme is linked and activated on site pages.
<!-- albersdown:theme-note:end -->
