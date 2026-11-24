# Outbreak probabilities for a spatial branching process

Computes the maximal fixed point of the multitype branching-process
equations in Wang and Valdano (2026), using the convention that rows of
`R` are destinations and columns are sources. Iteration starts at a
vector of ones. The Poisson calculation uses `-expm1(-crossprod(R, p))`.
The negative-binomial calculation implements Eq. (5) with
[`log1p()`](https://rdrr.io/r/base/Log.html) and
[`expm1()`](https://rdrr.io/r/base/Log.html) for numerical stability.

## Usage

``` r
outbreak_prob(
  R,
  tol = 1e-12,
  max_iter = 1e+05,
  offspring = "poisson",
  omega = 0
)
```

## Arguments

- R:

  A finite, nonnegative square reproduction matrix. `R[i, j]` is the
  expected number of infections in destination `i` caused by one
  infected individual from source `j`.

- tol:

  Positive convergence tolerance for both the iteration change and
  fixed-point residual.

- max_iter:

  Positive integer giving the maximum number of iterations.

- offspring:

  Offspring distribution. Either `"poisson"` (the default) or `"nb"` for
  independent negative-binomial counts on directed edges.

- omega:

  A finite nonnegative scalar overdispersion parameter. For the
  negative-binomial model, `size = 1 / omega` and
  `Var(X[i,j]) = R[i,j] * (1 + omega * R[i,j])`. The value zero recovers
  the Poisson model. Heterogeneous dispersion matrices are unsupported.

## Value

A named numeric vector, indexed by source community.

## Details

In subcritical and critical regimes, the theoretically correct maximal
fixed point is zero. The function returns it directly when the spectral
radius is at or below one (within a relative numerical tolerance),
avoiding arbitrarily slow convergence near the threshold. Above the
threshold, monotone iteration is supplemented by safeguarded Newton
steps when useful; the fixed-point residual is still checked before a
result is returned.

## Examples

``` r
R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
            dimnames = list(c("A", "B"), c("A", "B")))
outbreak_prob(R)
#>         A         B 
#> 0.4826325 0.2661426 
```
