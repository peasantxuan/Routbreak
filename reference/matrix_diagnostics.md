# Diagnose connected components and periodicity

Reports strongly connected components of the directed reproduction
graph. `cyclic_period = 1` denotes a primitive component; values greater
than one denote periodic (imprimitive) components. For an isolated
vertex without a self-loop, the cyclic period is undefined and reported
as `NA`.

## Usage

``` r
matrix_diagnostics(R)
```

## Arguments

- R:

  A finite, nonnegative square reproduction matrix. `R[i, j]` is the
  expected number of infections in destination `i` caused by one
  infected individual from source `j`.

## Value

A data frame with one row per community and its component number,
component size, cyclic period, primitive status, and matrix-level
irreducibility status.

## Details

Periodicity is separate from reducibility. In particular, the published
generation-limit continuation used by
[`r_ob()`](https://peasantxuan.github.io/Routbreak/reference/r_ob.md)
need not exist for an irreducible but periodic subcritical or critical
matrix. Routbreak reports that case rather than replacing the limit with
a cycle average.

## Examples

``` r
R <- matrix(c(1, 0.2, 0.1, 1), 2, byrow = TRUE,
            dimnames = list(c("A", "B"), c("A", "B")))
matrix_diagnostics(R)
#>   community component component_size cyclic_period primitive_component
#> 1         A         1              2             1                TRUE
#> 2         B         1              2             1                TRUE
#>   matrix_irreducible
#> 1               TRUE
#> 2               TRUE
```
