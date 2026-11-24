# Test whether a reproduction matrix is irreducible

A matrix is irreducible when its directed transmission graph is strongly
connected: every community can reach every other community through one
or more positive-weight edges.

## Usage

``` r
is_irreducible(R)
```

## Arguments

- R:

  A finite, nonnegative square reproduction matrix. `R[i, j]` is the
  expected number of infections in destination `i` caused by one
  infected individual from source `j`.

## Value

A single logical value.

## Examples

``` r
R <- matrix(c(1, 0.2, 0.1, 1), 2, byrow = TRUE,
            dimnames = list(c("A", "B"), c("A", "B")))
is_irreducible(R)
#> [1] TRUE
```
