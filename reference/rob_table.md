# Summarize outbreak risk by community

Calculates outbreak probabilities once and reuses the validated result
to report the principal community-level quantities. Rows of `R` are
destinations and columns are sources.

## Usage

``` r
rob_table(R, offspring = "poisson", omega = 0)
```

## Arguments

- R:

  A finite, nonnegative square reproduction matrix. `R[i, j]` is the
  expected number of infections in destination `i` caused by one
  infected individual from source `j`.

- offspring:

  Offspring distribution. Either `"poisson"` (the default) or `"nb"` for
  independent negative-binomial counts on directed edges.

- omega:

  A finite nonnegative scalar overdispersion parameter. For the
  negative-binomial model, `size = 1 / omega` and
  `Var(X[i,j]) = R[i,j] * (1 + omega * R[i,j])`. The value zero recovers
  the Poisson model. Heterogeneous dispersion matrices are unsupported.

## Value

A data frame with one row per source community and columns `community`,
`outbreak_probability`, `R_ob`, `R_local`, `R_total`, and `R_ref`.

## Examples

``` r
R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
            dimnames = list(c("A", "B"), c("A", "B")))
rob_table(R)
#>   community outbreak_probability     R_ob R_local R_total    R_ref
#> 1         A            0.4826325 1.365432     1.2     1.5 1.316228
#> 2         B            0.2661426 1.162687     0.8     1.0 1.316228
```
