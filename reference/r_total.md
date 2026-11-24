# Total reproduction ratios

Returns the expected total number of infections generated anywhere by an
infected individual from each source community.

## Usage

``` r
r_total(R)
```

## Arguments

- R:

  A finite, nonnegative square reproduction matrix. `R[i, j]` is the
  expected number of infections in destination `i` caused by one
  infected individual from source `j`.

## Value

A named numeric vector equal to `colSums(R)`.

## Examples

``` r
R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
            dimnames = list(c("A", "B"), c("A", "B")))
r_total(R)
#>   A   B 
#> 1.5 1.0 
```
