# Changelog

## Routbreak 0.3.0

- Added
  [`rob_table()`](https://peasantxuan.github.io/Routbreak/reference/rob_table.md)
  for a single aligned community-level summary and
  [`as_rmatrix()`](https://peasantxuan.github.io/Routbreak/reference/as_rmatrix.md)
  for validated destination-row/source-column long data.
- Added
  [`is_irreducible()`](https://peasantxuan.github.io/Routbreak/reference/is_irreducible.md)
  and
  [`matrix_diagnostics()`](https://peasantxuan.github.io/Routbreak/reference/matrix_diagnostics.md)
  to expose strongly connected components and cyclic periods. Periodic
  continuation cases are reported rather than assigned a cycle-averaged
  value.
- Replaced the dimension-dependent per-edge Gamma default with an
  explicit prior interface:
  [`infer_R()`](https://peasantxuan.github.io/Routbreak/reference/infer_R.md)
  and
  [`infer_R_sufficient()`](https://peasantxuan.github.io/Routbreak/reference/infer_R_sufficient.md)
  now require users to supply scientifically justified `prior_shape` and
  `prior_rate` values.
- Added
  [`infer_R()`](https://peasantxuan.github.io/Routbreak/reference/infer_R.md)
  for generation-resolved directed transmission data and
  [`infer_R_sufficient()`](https://peasantxuan.github.io/Routbreak/reference/infer_R_sufficient.md)
  for precomputed transmission counts and source exposures.
- Added explicit maximum-likelihood estimates, identification status,
  and proper shape-rate Gamma posterior moments for every transmission
  edge.
- Added
  [`estimate_rob()`](https://peasantxuan.github.io/Routbreak/reference/estimate_rob.md)
  to propagate posterior reproduction-matrix uncertainty to outbreak
  probabilities, outbreak reproduction ratios, and the reference
  reproduction ratio, with credible intervals and draw-failure
  diagnostics.
- Added the “Estimating outbreak reproduction ratios from early
  transmission data” vignette and deterministic and seeded inference
  tests.
- Improved near-threshold fixed-point convergence checks and prevented
  Gamma sampling underflow from creating artificial structural zeros.
- Documented sensitivity of sparse estimates to independent-entry priors
  and distinguished independent-entry Gamma–Poisson inference from the
  unimplemented Canadian structured Stan model `R = R_ref * M`.
- The inference likelihood remains Poisson-only. Negative-binomial
  zero-probability continuation and general reducible non-diagonal
  [`r_ob()`](https://peasantxuan.github.io/Routbreak/reference/r_ob.md)
  calculations remain explicitly unsupported rather than silently
  approximated.

## Routbreak 0.2.0

- Added independent-edge negative-binomial offspring distributions to
  [`outbreak_prob()`](https://peasantxuan.github.io/Routbreak/reference/outbreak_prob.md)
  and
  [`r_ob()`](https://peasantxuan.github.io/Routbreak/reference/r_ob.md)
  through the new `offspring` and `omega` arguments.
- Implemented the negative-binomial fixed point from Eq. (5) and the
  outbreak reproduction ratio from Eq. (S4) of Wang and Valdano (2026).
- Added deterministic, Poisson-limit, orientation, convergence, and
  seeded branching-process simulation tests.
- Added the “Overdispersed offspring distributions” vignette.

## Routbreak 0.1.0

- Initial mathematical implementation for multitype Poisson branching
  processes.
- Added outbreak, reference, local, and total reproduction ratios and
  matrix preparation helpers.
