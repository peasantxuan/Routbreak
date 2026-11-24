make_one_generation_data <- function() {
  data.frame(
    generation = rep(1L, 4L),
    source = c("A", "B", "A", "B"),
    destination = c("A", "A", "B", "B"),
    transmissions = c(12L, 4L, 3L, 16L),
    stringsAsFactors = FALSE
  )
}

test_that("generation-resolved inference has the published orientation", {
  fit <- infer_R(
    make_one_generation_data(),
    initial_cases = c(A = 10, B = 20),
    prior_shape = 2,
    prior_rate = 3
  )

  expected_S <- matrix(
    c(12, 4, 3, 16), nrow = 2, byrow = TRUE,
    dimnames = list(c("A", "B"), c("A", "B"))
  )
  expected_mle <- matrix(
    c(1.2, 0.2, 0.3, 0.8), nrow = 2, byrow = TRUE,
    dimnames = dimnames(expected_S)
  )

  expect_equal(fit$S, expected_S)
  expect_equal(fit$E, c(A = 10, B = 20))
  expect_equal(fit$mle, expected_mle)
  expect_equal(fit$alpha_post, expected_S + 2)
  expect_equal(fit$beta_post, matrix(
    c(13, 13, 23, 23), 2, 2,
    dimnames = dimnames(expected_S)
  ))
  expect_equal(fit$posterior_mean, fit$alpha_post / fit$beta_post)
  expect_equal(fit$posterior_variance,
               fit$alpha_post / fit$beta_post^2)
})

test_that("exposure uses the required previous-generation time shift", {
  data <- data.frame(
    generation = c(rep(1L, 4), rep(2L, 4)),
    source = rep(c("A", "B", "A", "B"), 2),
    destination = rep(c("A", "A", "B", "B"), 2),
    transmissions = c(1, 0, 1, 1, 1, 1, 0, 1),
    stringsAsFactors = FALSE
  )
  fit <- infer_R(data, c(A = 2, B = 1), prior_shape = 1, prior_rate = 1)

  expect_equal(
    fit$infections_by_generation,
    matrix(c(2, 1, 1, 2, 2, 1), nrow = 3, byrow = TRUE,
           dimnames = list(generation = c("0", "1", "2"),
                           community = c("A", "B")))
  )
  expect_equal(fit$E, c(A = 3, B = 3))
  expect_equal(
    fit$S,
    matrix(c(2, 1, 1, 2), 2, byrow = TRUE,
           dimnames = list(c("A", "B"), c("A", "B")))
  )
})

test_that("one-community inference and zero-exposure generation checks work", {
  one <- data.frame(
    generation = 1L, source = "A", destination = "A", transmissions = 3L
  )
  fit <- infer_R(one, c(A = 2), prior_shape = 1, prior_rate = 1)
  expect_equal(fit$S, matrix(3, 1, 1, dimnames = list("A", "A")))
  expect_equal(fit$E, c(A = 2))
  expect_equal(fit$mle, matrix(1.5, 1, 1,
                               dimnames = list("A", "A")))

  invalid <- data.frame(
    generation = 1L,
    source = "B",
    destination = "A",
    transmissions = 1L
  )
  expect_error(
    infer_R(invalid, c(A = 1, B = 0), prior_shape = 1, prior_rate = 1),
    "zero previous-generation exposure"
  )
})

test_that("zero exposure is explicit and inconsistent data are rejected", {
  S <- matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  fit <- infer_R_sufficient(
    S, c(A = 0, B = 5), prior_shape = 1, prior_rate = 1
  )

  expect_true(all(is.na(fit$mle[, "A"])))
  expect_equal(unname(fit$status[, "A"]), c("unidentified", "unidentified"))
  expect_true(all(is.finite(fit$posterior_mean)))

  S["B", "A"] <- 1
  expect_error(
    infer_R_sufficient(
      S, c(A = 0, B = 5), prior_shape = 1, prior_rate = 1
    ),
    "zero exposure"
  )

  bad <- make_one_generation_data()
  bad <- rbind(bad, bad[1, ])
  expect_error(
    infer_R(bad, c(A = 10, B = 20), prior_shape = 1, prior_rate = 1),
    "duplicate"
  )

  bad <- make_one_generation_data()
  bad$source[1] <- "C"
  expect_error(
    infer_R(bad, c(A = 10, B = 20), prior_shape = 1, prior_rate = 1),
    "absent"
  )
})

test_that("Gamma priors and structural zeros are validated", {
  ids <- c("A", "B")
  S <- matrix(c(2, 1, 0, 3), 2, byrow = TRUE,
              dimnames = list(ids, ids))
  prior_shape <- matrix(c(1, 2, 3, 4), 2, dimnames = list(ids, ids))
  prior_rate <- matrix(c(4, 3, 2, 1), 2, dimnames = list(ids, ids))
  zeros <- matrix(FALSE, 2, 2, dimnames = list(ids, ids))
  zeros[2, 1] <- TRUE
  fit <- infer_R_sufficient(
    S, c(A = 4, B = 6), prior_shape, prior_rate, zeros
  )

  expect_equal(fit$mle[2, 1], 0)
  expect_equal(fit$status[2, 1], "structural_zero")
  expect_equal(fit$posterior_mean[2, 1], 0)
  expect_true(is.na(fit$alpha_post[2, 1]))

  S[2, 1] <- 1
  expect_error(
    infer_R_sufficient(
      S, c(A = 4, B = 6), prior_shape = 1, prior_rate = 1,
      structural_zero = zeros
    ),
    "structural_zero"
  )
  expect_error(infer_R_sufficient(
    S, c(A = 4, B = 6), prior_shape = 0, prior_rate = 1
  ),
               "positive")
})

test_that("Gamma priors must be explicit and retain their intended scale", {
  ids <- c("A", "B")
  S <- matrix(0, 2, 2, dimnames = list(ids, ids))
  E <- c(A = 5, B = 5)
  expect_error(infer_R_sufficient(S, E), "Specify both")
  expect_error(
    infer_R_sufficient(S, E, prior_shape = 1),
    "Specify both"
  )

  # Equal edge weights can be elicited from a prior on each source total.
  total_shape <- 2
  total_rate <- 1.6
  edge_shape <- matrix(total_shape / length(ids), 2, 2,
                       dimnames = list(ids, ids))
  edge_rate <- matrix(total_rate, 2, 2, dimnames = list(ids, ids))
  fit <- infer_R_sufficient(
    S, E, prior_shape = edge_shape, prior_rate = edge_rate
  )
  expect_equal(
    colSums(fit$prior_shape / fit$prior_rate),
    c(A = total_shape / total_rate, B = total_shape / total_rate)
  )

  # A universal per-edge prior changes column-total scale with dimension.
  implied_total_mean <- function(n, shape = 0.001, rate = 0.001) {
    n * shape / rate
  }
  expect_equal(
    vapply(c(2, 10, 50, 110), implied_total_mean, numeric(1)),
    c(2, 10, 50, 110)
  )
})

test_that("posterior outbreak estimates are seeded and preserve identifiers", {
  fit <- infer_R(
    make_one_generation_data(), c(A = 10, B = 20),
    prior_shape = 1, prior_rate = 1
  )
  first <- estimate_rob(fit, draws = 30, seed = 0, return_draws = TRUE)
  second <- estimate_rob(fit, draws = 30, seed = 0, return_draws = TRUE)

  expect_equal(first$summary, second$summary)
  expect_equal(first$matrix_draws, second$matrix_draws)
  expect_equal(first$successful_draws, 30)
  expect_equal(first$failed_draws, 0)
  expect_equal(first$summary$community, c("A", "B"))
  expect_identical(dim(first$matrix_draws), c(30L, 2L, 2L))
  expect_true(all(first$matrix_draws >= 0))
  expect_true(all(first$summary$R_ob_lower <= first$summary$R_ob_median))
  expect_true(all(first$summary$R_ob_median <= first$summary$R_ob_upper))
})

test_that("seeded posterior sampling restores an existing RNG state", {
  fit <- infer_R(
    make_one_generation_data(), c(A = 10, B = 20),
    prior_shape = 1, prior_rate = 1
  )
  set.seed(20261124)
  state_before <- .Random.seed

  first <- estimate_rob(fit, draws = 20, seed = 7, return_draws = TRUE)
  expect_identical(.Random.seed, state_before)
  second <- estimate_rob(fit, draws = 20, seed = 7, return_draws = TRUE)

  expect_identical(.Random.seed, state_before)
  expect_equal(first$matrix_draws, second$matrix_draws)
})

test_that("seeded posterior sampling does not leave a new RNG state", {
  fit <- infer_R(
    make_one_generation_data(), c(A = 10, B = 20),
    prior_shape = 1, prior_rate = 1
  )
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    saved_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  estimate_rob(fit, draws = 20, seed = 8)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("seeded posterior sampling restores RNG state after errors", {
  ids <- c("A", "B")
  S <- matrix(c(4, 0, 2, 4), 2, byrow = TRUE,
              dimnames = list(ids, ids))
  zeros <- matrix(FALSE, 2, 2, dimnames = list(ids, ids))
  zeros[1, 2] <- TRUE
  fit <- infer_R_sufficient(
    S, c(A = 2, B = 2), prior_shape = 1, prior_rate = 1,
    structural_zero = zeros
  )
  set.seed(20261125)
  state_before <- .Random.seed

  expect_error(
    estimate_rob(
      fit, draws = 5, seed = 9, max_fail_fraction = 0.9
    ),
    "posterior draws failed"
  )
  expect_identical(.Random.seed, state_before)
})

test_that("credible intervals are empirical posterior quantiles", {
  S <- matrix(12, 1, 1, dimnames = list("A", "A"))
  fit <- infer_R_sufficient(S, c(A = 10), prior_shape = 2, prior_rate = 3)
  level <- 0.8
  estimate <- estimate_rob(
    fit, draws = 200, seed = 42, level = level, return_draws = TRUE
  )
  probabilities <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)

  expect_equal(
    as.numeric(estimate$summary[1, c(
      "R_ob_lower", "R_ob_median", "R_ob_upper"
    )]),
    unname(quantile(estimate$R_ob_draws[, 1], probabilities))
  )
  expect_equal(
    as.numeric(estimate$summary[1, c(
      "outbreak_probability_lower",
      "outbreak_probability_median",
      "outbreak_probability_upper"
    )]),
    unname(quantile(
      estimate$outbreak_probability_draws[, 1], probabilities
    ))
  )
  expect_equal(
    as.numeric(estimate$summary[1, c(
      "R_ref_lower", "R_ref_median", "R_ref_upper"
    )]),
    unname(quantile(estimate$R_ref_draws, probabilities))
  )
})

test_that("Gamma underflow does not create artificial structural zeros", {
  ids <- c("A", "B")
  S <- matrix(0, 2, 2, dimnames = list(ids, ids))
  fit <- infer_R_sufficient(
    S, c(A = 1, B = 1), prior_shape = 0.001, prior_rate = 0.001
  )
  estimate <- estimate_rob(
    fit, draws = 20, seed = 1, return_draws = TRUE
  )

  expect_equal(estimate$successful_draws, 20)
  expect_equal(estimate$failed_draws, 0)
  expect_gt(estimate$diagnostics$gamma_underflow_replacements, 0)
  expect_true(all(estimate$matrix_draws > 0))
})

test_that("declared structural zeros remain exact in posterior draws", {
  ids <- c("A", "B")
  S <- matrix(c(8, 0, 0, 7), 2, byrow = TRUE,
              dimnames = list(ids, ids))
  structural_zero <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2, byrow = TRUE,
                            dimnames = list(ids, ids))
  fit <- infer_R_sufficient(
    S, c(A = 5, B = 5), prior_shape = 1, prior_rate = 1,
    structural_zero = structural_zero
  )
  estimate <- estimate_rob(
    fit, draws = 20, seed = 9, return_draws = TRUE
  )

  expect_true(all(estimate$matrix_draws[, 1, 2] == 0))
  expect_true(all(estimate$matrix_draws[, 2, 1] == 0))
  expect_true(all(estimate$matrix_draws[, 1, 1] > 0))
  expect_true(all(estimate$matrix_draws[, 2, 2] > 0))
  expect_equal(estimate$failed_draws, 0)
})

test_that("unsupported posterior matrix structures are not discarded", {
  ids <- c("A", "B")
  S <- matrix(c(4, 0, 2, 4), 2, byrow = TRUE,
              dimnames = list(ids, ids))
  zeros <- matrix(FALSE, 2, 2, dimnames = list(ids, ids))
  zeros[1, 2] <- TRUE
  fit <- infer_R_sufficient(
    S, c(A = 2, B = 2), prior_shape = 1, prior_rate = 1,
    structural_zero = zeros
  )

  expect_error(
    estimate_rob(fit, draws = 5, seed = 11, max_fail_fraction = 0.9),
    "posterior draws failed"
  )
})

test_that("a seeded high-exposure simulation recovers matrix orientation", {
  R_true <- matrix(
    c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
    dimnames = list(c("A", "B"), c("A", "B"))
  )
  initial <- c(A = 5000, B = 5000)
  set.seed(814)
  Y <- matrix(
    rpois(4, lambda = as.vector(sweep(R_true, 2, initial, `*`))),
    2, 2, dimnames = dimnames(R_true)
  )
  data <- expand.grid(
    destination = rownames(R_true), source = colnames(R_true),
    stringsAsFactors = FALSE
  )
  data$generation <- 1L
  data$transmissions <- mapply(
    function(destination, source) Y[destination, source],
    data$destination, data$source
  )
  data <- data[, c("generation", "source", "destination", "transmissions")]
  fit <- infer_R(data, initial, prior_shape = 1, prior_rate = 1)

  expect_equal(fit$mle, R_true, tolerance = 0.04)
})

test_that("inference controls reject invalid values", {
  fit <- infer_R(
    make_one_generation_data(), c(A = 10, B = 20),
    prior_shape = 1, prior_rate = 1
  )
  expect_error(estimate_rob(fit, draws = 0), "positive integer")
  expect_error(estimate_rob(fit, seed = -1), "nonnegative integer")
  expect_error(estimate_rob(fit, level = 1), "strictly between")
  expect_error(estimate_rob(fit, max_fail_fraction = 1), "max_fail_fraction")
})
