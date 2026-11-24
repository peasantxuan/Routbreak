# Package index

## Core outbreak metrics

Introduction-specific early-epidemic risk.

- [`outbreak_prob()`](https://peasantxuan.github.io/Routbreak/reference/outbreak_prob.md)
  : Outbreak probabilities for a spatial branching process
- [`r_ob()`](https://peasantxuan.github.io/Routbreak/reference/r_ob.md)
  : Outbreak reproduction ratio
- [`r_ref()`](https://peasantxuan.github.io/Routbreak/reference/r_ref.md)
  : Reference reproduction ratio
- [`r_local()`](https://peasantxuan.github.io/Routbreak/reference/r_local.md)
  : Local reproduction ratios
- [`r_total()`](https://peasantxuan.github.io/Routbreak/reference/r_total.md)
  : Total reproduction ratios
- [`rob_table()`](https://peasantxuan.github.io/Routbreak/reference/rob_table.md)
  : Summarize outbreak risk by community

## Matrix utilities

Prepare, scale, and diagnose reproduction matrices.

- [`scale_to_rref()`](https://peasantxuan.github.io/Routbreak/reference/scale_to_rref.md)
  : Scale a reproduction matrix to a target reference ratio
- [`colocation_to_r()`](https://peasantxuan.github.io/Routbreak/reference/colocation_to_r.md)
  : Construct a reproduction matrix from colocation data
- [`as_rmatrix()`](https://peasantxuan.github.io/Routbreak/reference/as_rmatrix.md)
  : Convert directed long-format data to a reproduction matrix
- [`is_irreducible()`](https://peasantxuan.github.io/Routbreak/reference/is_irreducible.md)
  : Test whether a reproduction matrix is irreducible
- [`matrix_diagnostics()`](https://peasantxuan.github.io/Routbreak/reference/matrix_diagnostics.md)
  : Diagnose connected components and periodicity

## Statistical inference

Estimate reproduction matrices and propagate posterior uncertainty.

- [`infer_R()`](https://peasantxuan.github.io/Routbreak/reference/infer_R.md)
  : Infer a reproduction matrix from generation-resolved transmission
  data
- [`infer_R_sufficient()`](https://peasantxuan.github.io/Routbreak/reference/infer_R_sufficient.md)
  : Infer a reproduction matrix from sufficient statistics
- [`estimate_rob()`](https://peasantxuan.github.io/Routbreak/reference/estimate_rob.md)
  : Propagate reproduction-matrix uncertainty to outbreak risk
