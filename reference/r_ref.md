# Reference reproduction ratio

Reference reproduction ratio

## Usage

``` r
r_ref(R)
```

## Arguments

- R:

  A finite, nonnegative square reproduction matrix. `R[i, j]` is the
  expected number of infections in destination `i` caused by one
  infected individual from source `j`.

## Value

A single numeric value equal to the spectral radius of `R`.

## Examples

``` r
R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE)
r_ref(R)
#> [1] 1.316228
```
