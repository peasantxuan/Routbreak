# Scale a reproduction matrix to a target reference ratio

Scale a reproduction matrix to a target reference ratio

## Usage

``` r
scale_to_rref(R, target)
```

## Arguments

- R:

  A finite, nonnegative square reproduction matrix. `R[i, j]` is the
  expected number of infections in destination `i` caused by one
  infected individual from source `j`.

- target:

  A finite nonnegative target spectral radius.

## Value

A scaled reproduction matrix with the original dimnames.

## Examples

``` r
R <- matrix(c(1, 0.2, 0.1, 0.8), 2, byrow = TRUE)
scaled <- scale_to_rref(R, target = 1.5)
r_ref(scaled)
#> [1] 1.5
```
