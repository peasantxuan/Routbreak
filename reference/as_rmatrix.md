# Convert directed long-format data to a reproduction matrix

Creates a square matrix whose rows are destinations and columns are
sources. An explicit zero-valued row creates a zero matrix entry. An
absent source-destination pair remains an unobserved edge and causes an
error by default; it is converted to zero only when `missing = "zero"`
is requested.

## Usage

``` r
as_rmatrix(
  data,
  source,
  destination,
  value,
  communities = NULL,
  missing = c("error", "zero")
)
```

## Arguments

- data:

  A data frame containing directed edge data.

- source, destination, value:

  Single strings naming the source, destination, and nonnegative numeric
  value columns.

- communities:

  Optional character vector defining the complete set and order of
  community identifiers. When omitted, order is determined by the first
  appearance among destinations and then sources.

- missing:

  How absent directed edges are handled. The default, `"error"`, rejects
  incomplete edge data. Use `"zero"` only when absence is scientifically
  intended to mean a zero matrix entry.

## Value

A numeric reproduction matrix with destination rows, source columns, and
named `destination` and `source` dimname axes.

## Examples

``` r
edges <- data.frame(
  source = c("A", "A", "B", "B"),
  destination = c("A", "B", "A", "B"),
  value = c(1.2, 0.3, 0.2, 0.8)
)
as_rmatrix(edges, "source", "destination", "value")
#>            source
#> destination   A   B
#>           A 1.2 0.2
#>           B 0.3 0.8
```
