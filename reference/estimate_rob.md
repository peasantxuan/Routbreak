# Propagate reproduction-matrix uncertainty to outbreak risk

Draws complete reproduction matrices from a fitted Gamma posterior and
applies
[`outbreak_prob()`](https://peasantxuan.github.io/Routbreak/reference/outbreak_prob.md),
[`r_ob()`](https://peasantxuan.github.io/Routbreak/reference/r_ob.md),
and optionally
[`r_ref()`](https://peasantxuan.github.io/Routbreak/reference/r_ref.md)
to every draw. The transmission likelihood is Poisson; negative-binomial
transmission inference is not implied by the offspring-distribution
option in
[`outbreak_prob()`](https://peasantxuan.github.io/Routbreak/reference/outbreak_prob.md).
Intervals are posterior credible intervals, not frequentist bootstrap
intervals.

## Usage

``` r
estimate_rob(
  fit,
  draws = 1000,
  seed = NULL,
  level = 0.95,
  return_draws = FALSE,
  include_rref = TRUE,
  max_fail_fraction = 0.05,
  tol = 1e-12,
  max_iter = 1e+05
)
```

## Arguments

- fit:

  An object returned by
  [`infer_R()`](https://peasantxuan.github.io/Routbreak/reference/infer_R.md)
  or
  [`infer_R_sufficient()`](https://peasantxuan.github.io/Routbreak/reference/infer_R_sufficient.md).

- draws:

  Positive integer number of posterior matrices to sample.

- seed:

  Optional finite integer random seed.

- level:

  Credible interval level strictly between zero and one.

- return_draws:

  If `TRUE`, include posterior matrix, outbreak-ratio, and
  outbreak-probability draws in the returned object.

- include_rref:

  If `TRUE`, calculate the reference reproduction ratio for every
  successful draw.

- max_fail_fraction:

  Maximum tolerated fraction of failed draws.

- tol, max_iter:

  Numerical controls passed to
  [`outbreak_prob()`](https://peasantxuan.github.io/Routbreak/reference/outbreak_prob.md).

## Value

An object of class `routbreak_rob_estimate`. Its `summary` data frame
contains posterior medians and lower/upper credible limits for `R_ob`,
outbreak probability, and (when requested) `R_ref`. `requested_draws`,
`successful_draws`, `failed_draws`, `failures`, and `diagnostics` report
the failure fraction, reasons, numerical controls, and Gamma-underflow
replacements. If `return_draws = TRUE`, successful matrix, `R_ob`,
outbreak-probability, and optional `R_ref` draws are returned in arrays
or matrices preserving community identifiers.

## Details

Draws that fail numerical or theoretical checks are recorded. Any
failures produce a warning, and the function errors when their fraction
exceeds `max_fail_fraction`. R's Gamma generator can underflow to exact
zero for very small shape parameters even though a proper Gamma
distribution has strictly positive support. Such values are replaced by
the smallest positive normalized floating-point value and counted in the
diagnostics; this prevents numerical underflow from creating artificial
structural zeros.

Posterior draws are independent Gamma samples, not an MCMC chain, so
MCMC convergence diagnostics do not apply. Reported quantiles still have
Monte Carlo error; users should increase `draws` and check stability for
final analyses. If some draws fail within the allowed fraction,
summaries are conditional on the successful draws and must be
interpreted with the accompanying failure diagnostics.

## Examples

``` r
S <- matrix(12, 1, 1, dimnames = list("A", "A"))
fit <- infer_R_sufficient(
  S, E = c(A = 10), prior_shape = 1, prior_rate = 1
)
estimate <- estimate_rob(fit, draws = 20, seed = 1)
estimate$summary
#>   community R_ob_median R_ob_lower R_ob_upper outbreak_probability_median
#> A         A    1.216999  0.6383935   1.593597                   0.3297678
#>   outbreak_probability_lower outbreak_probability_upper R_ref_median
#> A                          0                  0.6384559     1.216999
#>   R_ref_lower R_ref_upper
#> A   0.6383935    1.593597
estimate$diagnostics
#> $tolerance
#> [1] 1e-12
#> 
#> $max_iterations
#> [1] 100000
#> 
#> $failure_fraction
#> [1] 0
#> 
#> $gamma_underflow_replacements
#> [1] 0
#> 
#> $convergence
#> [1] "Every successful draw satisfied the fixed-point residual check implemented by outbreak_prob()."
#> 
```
