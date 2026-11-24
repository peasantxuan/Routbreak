named_matrix <- function(x, ids = paste0("c", seq_len(nrow(x)))) {
  dimnames(x) <- list(ids, ids)
  x
}

test_that("single population benchmark is correct", {
  R <- named_matrix(matrix(1.5, 1, 1), "one")
  expect_equal(outbreak_prob(R), c(one = 0.5828116438657147), tolerance = 1e-8)
  expect_equal(r_ob(R), c(one = 1.5), tolerance = 1e-12)
  expect_equal(r_ref(R), 1.5, tolerance = 1e-12)
  expect_equal(r_local(R), c(one = 1.5), tolerance = 1e-12)
  expect_equal(r_total(R), c(one = 1.5), tolerance = 1e-12)
})

test_that("independent communities are treated separately", {
  R <- named_matrix(matrix(c(0.5, 0, 0, 2), 2, byrow = TRUE))
  expect_equal(outbreak_prob(R), c(c1 = 0, c2 = 0.79681213002002),
               tolerance = 1e-8)
  expect_equal(r_ob(R), c(c1 = 0.5, c2 = 2), tolerance = 1e-12)
  expect_equal(r_ref(R), 2, tolerance = 1e-12)
})

test_that("homogeneous column totals give homogeneous outbreak ratios", {
  R <- named_matrix(matrix(c(0.8, 0.3, 0.4, 0.9), 2, byrow = TRUE))
  expect_equal(r_ref(R), 1.2, tolerance = 1e-12)
  expect_equal(outbreak_prob(R), c(c1 = 0.31369833104122, c2 = 0.31369833104122),
               tolerance = 1e-8)
  expect_equal(r_ob(R), c(c1 = 1.2, c2 = 1.2), tolerance = 1e-8)
  expect_equal(r_total(R), c(c1 = 1.2, c2 = 1.2), tolerance = 1e-12)
})

test_that("asymmetric benchmark fixes matrix orientation", {
  R <- named_matrix(matrix(c(1.3, 0.05, 0.5, 0.8), 2, byrow = TRUE))
  expect_equal(outbreak_prob(R),
               c(c1 = 0.50980927634099, c2 = 0.10041734481623),
               tolerance = 1e-8)
  expect_equal(r_ob(R),
               c(c1 = 1.39848520758287, c2 = 1.05384522826905),
               tolerance = 1e-8)
})

test_that("irreducible subcritical continuation equals spectral radius", {
  R <- named_matrix(matrix(c(0.4, 0.2, 0.05, 0.6), 2, byrow = TRUE))
  expected <- 0.6414213562373095
  expect_equal(outbreak_prob(R), c(c1 = 0, c2 = 0), tolerance = 1e-12)
  expect_equal(r_ref(R), expected, tolerance = 1e-12)
  expect_equal(r_ob(R), c(c1 = expected, c2 = expected), tolerance = 1e-12)
})

test_that("critical irreducible matrices use the theoretical limit", {
  R <- named_matrix(matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE))
  expect_equal(outbreak_prob(R), c(c1 = 0, c2 = 0), tolerance = 0)
  expect_equal(r_ob(R), c(c1 = 1, c2 = 1), tolerance = 1e-12)
})

test_that("near-critical probabilities are resolved without false failure", {
  subcritical <- matrix(0.99999, 1, 1,
                        dimnames = list("one", "one"))
  expect_equal(
    outbreak_prob(subcritical, max_iter = 1),
    c(one = 0),
    tolerance = 0
  )

  supercritical <- matrix(1.00001, 1, 1,
                          dimnames = list("one", "one"))
  independent_root <- uniroot(
    function(p) -expm1(-1.00001 * p) - p,
    c(.Machine$double.eps, 0.1),
    tol = .Machine$double.eps
  )$root
  expect_lt(
    abs(unname(outbreak_prob(supercritical)) - independent_root),
    1e-12
  )
})

test_that("colocation construction weights rows and reaches target", {
  ids <- c("north", "south")
  C <- matrix(c(0.1, 0.7, 0.2, 0.4), 2, byrow = TRUE,
              dimnames = list(ids, ids))
  population <- c(north = 100, south = 1000)
  target <- 2.5
  R <- colocation_to_r(C, population, target)
  A <- C * population
  expected <- target * A / r_ref(A)
  expect_equal(R, expected, tolerance = 1e-14)
  expect_equal(R[1, ] / C[1, ],
               setNames(rep(100 * target / r_ref(A), 2), ids),
               tolerance = 1e-14)
  expect_equal(R[2, ] / C[2, ],
               setNames(rep(1000 * target / r_ref(A), 2), ids),
               tolerance = 1e-14)
  expect_equal(r_ref(R), target, tolerance = 1e-12)
  expect_identical(dimnames(R), dimnames(C))
})

test_that("scaling preserves dimnames and reaches the target", {
  R <- named_matrix(matrix(c(1, 2, 3, 4), 2))
  scaled <- scale_to_rref(R, 3)
  expect_equal(r_ref(scaled), 3, tolerance = 1e-12)
  expect_identical(dimnames(scaled), dimnames(R))
  expect_equal(scale_to_rref(R, 0), R * 0, tolerance = 0)
})

test_that("invalid matrices and parameters fail clearly", {
  expect_error(outbreak_prob(data.frame(x = 1)), "numeric matrix")
  expect_error(outbreak_prob(matrix(1:6, 2, 3)), "square")
  expect_error(outbreak_prob(matrix(c(1, -0.1, 0, 1), 2)), "nonnegative")
  expect_error(outbreak_prob(matrix(c(1, NA, 0, 1), 2)), "finite")
  expect_error(outbreak_prob(matrix(1.5), tol = 0), "positive")
  expect_error(outbreak_prob(matrix(1.5), tol = "small"), "positive")
  expect_error(outbreak_prob(matrix(1.5), max_iter = 0), "positive integer")
  expect_error(outbreak_prob(matrix(1.5), max_iter = Inf), "positive integer")
  expect_error(outbreak_prob(matrix(1.5), max_iter = 1), "did not converge")
})

test_that("identifier mismatches and invalid populations fail", {
  C <- matrix(c(1, 0.2, 0.3, 1), 2,
              dimnames = list(c("a", "b"), c("a", "b")))
  expect_error(colocation_to_r(C, c(b = 10, a = 20), 2),
               "same order")
  expect_error(colocation_to_r(unname(C), c(a = 10, b = 20), 2),
               "identifiers")
  expect_error(colocation_to_r(C, c(a = 10, a = 20), 2), "unique")
  expect_error(colocation_to_r(C, c(a = 10, b = NA), 2), "finite positive")

  bad_names <- C
  colnames(bad_names) <- c("b", "a")
  expect_error(r_ref(bad_names), "same order")
})

test_that("zero spectral radius and unsupported reducibility are explicit", {
  zero <- named_matrix(matrix(0, 2, 2))
  expect_equal(outbreak_prob(zero), c(c1 = 0, c2 = 0), tolerance = 0)
  expect_equal(r_ob(zero), c(c1 = 0, c2 = 0), tolerance = 0)
  expect_error(scale_to_rref(zero, 1), "zero-spectral-radius")

  reducible <- named_matrix(matrix(c(1.2, 0.4, 0, 0.8), 2, byrow = TRUE))
  expect_error(r_ob(reducible), "reducible non-diagonal")
})

test_that("supplied probabilities are validated", {
  R <- named_matrix(matrix(c(0.8, 0.3, 0.4, 0.9), 2, byrow = TRUE))
  p <- outbreak_prob(R)
  expect_equal(r_ob(R, p), r_ob(R), tolerance = 1e-12)
  expect_error(r_ob(R, c(c2 = p[2], c1 = p[1])), "Names")
  expect_error(r_ob(R, c(c1 = 0, c2 = 0)), "maximal")
  expect_error(r_ob(R, c(c1 = 1.1, c2 = 0.2)), "between zero and one")
  expect_error(r_ob(R, p + 1e-9), "validation tolerance")
})

test_that("rob_table reports aligned community quantities and reuses p", {
  R <- named_matrix(
    matrix(c(1.3, 0.05, 0.5, 0.8), 2, byrow = TRUE),
    c("urban", "rural")
  )
  original_outbreak_prob <- outbreak_prob
  calls <- 0L
  testthat::local_mocked_bindings(
    outbreak_prob = function(...) {
      calls <<- calls + 1L
      original_outbreak_prob(...)
    },
    .package = "Routbreak"
  )

  result <- rob_table(R)
  expect_equal(calls, 1L)
  expect_identical(result$community, c("urban", "rural"))
  expect_equal(result$outbreak_probability, unname(original_outbreak_prob(R)))
  expect_equal(result$R_ob, unname(r_ob(R)))
  expect_equal(result$R_local, unname(diag(R)))
  expect_equal(result$R_total, unname(colSums(R)))
  expect_equal(result$R_ref, rep(r_ref(R), 2))
})

test_that("as_rmatrix preserves directed destination-row orientation", {
  edges <- data.frame(
    from = c("A", "B", "A", "B"),
    to = c("B", "A", "A", "B"),
    mean = c(0.3, 0.2, 1.1, 0.8)
  )
  R <- as_rmatrix(
    edges, source = "from", destination = "to", value = "mean",
    communities = c("A", "B")
  )
  expect_equal(
    R,
    matrix(c(1.1, 0.2, 0.3, 0.8), 2, byrow = TRUE,
           dimnames = list(
             destination = c("A", "B"), source = c("A", "B")
           ))
  )
  expect_identical(names(dimnames(R)), c("destination", "source"))
  expect_equal(R["B", "A"], 0.3)
})

test_that("as_rmatrix distinguishes absent edges from explicit zeros", {
  incomplete <- data.frame(
    source = c("A", "B", "A"),
    destination = c("A", "A", "B"),
    value = c(1, 0.2, 0)
  )
  expect_error(
    as_rmatrix(incomplete, "source", "destination", "value",
               communities = c("A", "B")),
    "incomplete"
  )
  completed <- as_rmatrix(
    incomplete, "source", "destination", "value",
    communities = c("A", "B"), missing = "zero"
  )
  expect_equal(completed["B", "A"], 0)
  expect_equal(completed["B", "B"], 0)

  duplicate <- rbind(incomplete, incomplete[1, ])
  expect_error(
    as_rmatrix(duplicate, "source", "destination", "value",
               communities = c("A", "B"), missing = "zero"),
    "duplicate"
  )
  expect_error(
    as_rmatrix(incomplete, "source", "destination", "value",
               communities = "A", missing = "zero"),
    "unknown"
  )
  incomplete$value[1] <- -1
  expect_error(
    as_rmatrix(incomplete, "source", "destination", "value",
               communities = c("A", "B"), missing = "zero"),
    "nonnegative"
  )
})

test_that("matrix diagnostics report components and cyclic period", {
  primitive <- named_matrix(matrix(c(0.4, 0.2, 0.1, 0.5), 2, byrow = TRUE))
  expect_true(is_irreducible(primitive))
  primitive_diagnostics <- matrix_diagnostics(primitive)
  expect_true(all(primitive_diagnostics$matrix_irreducible))
  expect_equal(unique(primitive_diagnostics$cyclic_period), 1L)
  expect_true(all(primitive_diagnostics$primitive_component))

  periodic <- named_matrix(matrix(c(0, 0.4, 0.6, 0), 2, byrow = TRUE))
  expect_true(is_irreducible(periodic))
  periodic_diagnostics <- matrix_diagnostics(periodic)
  expect_equal(unique(periodic_diagnostics$cyclic_period), 2L)
  expect_false(any(periodic_diagnostics$primitive_component))
  expect_error(r_ob(periodic), "periodic matrix")

  p_t <- c(1, 1)
  continuation_ratio <- matrix(NA_real_, 30, 2)
  for (generation in seq_len(30)) {
    continuation_ratio[generation, ] <-
      as.vector(crossprod(periodic, p_t)) / p_t
    p_t <- -expm1(-as.vector(crossprod(periodic, p_t)))
  }
  expect_equal(continuation_ratio[29, ],
               rev(continuation_ratio[30, ]), tolerance = 1e-7)
  expect_gt(max(abs(continuation_ratio[30, ] - r_ref(periodic))), 0.02)

  supercritical_periodic <- periodic * 3
  expect_gt(r_ref(supercritical_periodic), 1)
  expect_true(all(outbreak_prob(supercritical_periodic) > 0))
  expect_true(all(is.finite(r_ob(supercritical_periodic))))

  reducible <- named_matrix(matrix(c(0.5, 0, 0.2, 0.4), 2, byrow = TRUE))
  expect_false(is_irreducible(reducible))
  reducible_diagnostics <- matrix_diagnostics(reducible)
  expect_equal(length(unique(reducible_diagnostics$component)), 2L)
  expect_identical(reducible_diagnostics$community, c("c1", "c2"))
})

test_that("v0.1 positional arguments remain compatible", {
  R <- named_matrix(matrix(c(0.8, 0.3, 0.4, 0.9), 2, byrow = TRUE))
  expect_equal(
    outbreak_prob(R, 1e-10, 10000),
    outbreak_prob(R, tol = 1e-10, max_iter = 10000)
  )
  p <- outbreak_prob(R)
  expect_equal(r_ob(R, p), r_ob(R = R, p = p))
  expect_identical(names(formals(outbreak_prob))[1:3],
                   c("R", "tol", "max_iter"))
  expect_identical(names(formals(r_ob))[1:2], c("R", "p"))
})
