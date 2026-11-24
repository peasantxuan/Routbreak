# Infer a reproduction matrix from sufficient statistics

Fits the independent Poisson transmission model in Eqs. (27)–(30) of
Wang and Valdano (2026). `S[i, j]` is the total number of observed
transmissions from source `j` to destination `i`, and `E[j]` is the
total number of infectious source-`j` individuals exposed across the
preceding generations. Rows are destinations and columns are sources.

## Usage

``` r
infer_R_sufficient(
  S,
  E,
  prior_shape = NULL,
  prior_rate = NULL,
  structural_zero = NULL
)
```

## Arguments

- S:

  A finite nonnegative square numeric matrix of integer transmission
  counts. Row and column identifiers must be unique and identical in the
  same order.

- E:

  A named finite nonnegative numeric vector of integer source exposures,
  in the same order as the columns of `S`.

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

An object of class `routbreak_inference` containing sufficient
statistics `S` and `E`; square destination-by-source matrices `mle`,
`status`, `prior_shape`, `prior_rate`, `alpha_post`, `beta_post`,
`posterior_mean`, and `posterior_variance`; the logical matrix
`structural_zero`; and model metadata. The `mle` is `NA` for
unidentified zero-exposure entries, while their proper-prior posterior
moments remain defined. Declared structural zeros have zero MLE and
posterior moments.

## Details

The maximum-likelihood estimate is `S[i, j] / E[j]` when `E[j] > 0`.
Entries with zero exposure are reported as unidentified rather than set
to zero. Under a proper Gamma prior in shape-rate form, posterior shape
is `prior_shape + S` and posterior rate is `prior_rate + E[j]`.

Prior parameters are required explicitly because a fixed independent
prior on every edge implies a different prior for total reproduction as
the number of communities changes. For example, if every edge has a
common `Gamma(a, b)` prior, the prior mean of each source-column total
is `n * a / b`. Users should elicit priors for the matrix scale and
structure relevant to their study. This is not the paper's improper
flat-prior construction. The posterior mean is
`(prior_shape + S) / (prior_rate + E)`, whereas the MLE remains `S / E`.

## Examples

``` r
ids <- c("A", "B")
S <- matrix(c(12, 4, 3, 16), 2, byrow = TRUE,
            dimnames = list(ids, ids))
fit <- infer_R_sufficient(
  S, E = c(A = 10, B = 20), prior_shape = 1, prior_rate = 1
)
fit$mle
#>     A   B
#> A 1.2 0.2
#> B 0.3 0.8
fit$posterior_mean
#>           A         B
#> A 1.1818182 0.2380952
#> B 0.3636364 0.8095238
```
