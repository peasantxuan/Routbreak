# Construct a reproduction matrix from colocation data

Constructs `A[i, j] = C[i, j] * population[i]`, then scales `A` so its
spectral radius is `target_rref`. Community identifiers must be
explicit, unique, and identical across matrix rows, matrix columns, and
the population vector; inputs are never silently reordered.

## Usage

``` r
colocation_to_r(C, population, target_rref)
```

## Arguments

- C:

  A finite, nonnegative square numeric colocation matrix with matching
  row and column identifiers.

- population:

  A named finite positive numeric vector of destination population
  sizes.

- target_rref:

  A finite positive target reference reproduction ratio.

## Value

A reproduction matrix with the same dimnames as `C`.

## Examples

``` r
C <- matrix(c(0.8, 0.1, 0.2, 0.7), 2, byrow = TRUE,
            dimnames = list(c("A", "B"), c("A", "B")))
population <- c(A = 1000, B = 2000)
R <- colocation_to_r(C, population, target_rref = 1.5)
r_ref(R)
#> [1] 1.5
```
