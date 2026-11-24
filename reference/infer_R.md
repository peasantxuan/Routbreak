# Infer a reproduction matrix from generation-resolved transmission data

Converts generation-resolved directed transmission counts into the
sufficient statistics used by
[`infer_R_sufficient()`](https://peasantxuan.github.io/Routbreak/reference/infer_R_sufficient.md).
For generation `t`, `Y[t, i, j]` records infections in destination `i`
caused by source `j`. Source exposure at generation `t` is the number
infected in source `j` at generation `t - 1`; the initially infected
generation is therefore included in exposure for the first observed
transmission generation.

## Usage

``` r
infer_R(
  data,
  initial_cases,
  prior_shape = NULL,
  prior_rate = NULL,
  structural_zero = NULL
)
```

## Arguments

- data:

  A data frame with columns `generation`, `source`, `destination`, and
  `transmissions`. Generations are positive integers beginning with the
  first observed offspring generation.

- initial_cases:

  A named vector of finite nonnegative integer counts for generation
  zero. Its names define the canonical community order.

- prior_shape:

  A required finite positive scalar or matrix of Gamma prior shape
  parameters.

- prior_rate:

  A required finite positive scalar or matrix of Gamma prior rate
  parameters.

- structural_zero:

  Optional logical matrix identifying edges fixed at zero by design.
  Positive observations on those edges are rejected. A zero observed
  count is otherwise treated as a sampling zero, not a structural zero.

## Value

An object of class `routbreak_inference`. In addition to the fields
returned by
[`infer_R_sufficient()`](https://peasantxuan.github.io/Routbreak/reference/infer_R_sufficient.md),
it contains the reconstructed array `Y` with dimensions generation by
destination by source, named `initial_cases`, and the
generation-by-community matrix `infections_by_generation` used to
calculate previous-generation exposures.

## Details

Missing source-destination rows are interpreted as zero transmissions.
Duplicate generation-source-destination rows are rejected rather than
silently aggregated.

## Examples

``` r
transmissions <- data.frame(
  generation = rep(1L, 4),
  source = c("A", "B", "A", "B"),
  destination = c("A", "A", "B", "B"),
  transmissions = c(12L, 4L, 3L, 16L)
)
fit <- infer_R(
  transmissions, initial_cases = c(A = 10, B = 20),
  prior_shape = 1, prior_rate = 1
)
fit$S
#>    A  B
#> A 12  4
#> B  3 16
fit$E
#>  A  B 
#> 10 20 
fit$mle
#>     A   B
#> A 1.2 0.2
#> B 0.3 0.8
```
