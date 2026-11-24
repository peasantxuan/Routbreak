# Local reproduction ratios

Local reproduction ratios

## Usage

``` r
r_local(R)
```

## Arguments

- R:

  A finite, nonnegative square reproduction matrix. `R[i, j]` is the
  expected number of infections in destination `i` caused by one
  infected individual from source `j`.

## Value

A named numeric vector equal to `diag(R)`.

## Examples

``` r
R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
            dimnames = list(c("A", "B"), c("A", "B")))
r_local(R)
#>   A   B 
#> 1.2 0.8 
```
