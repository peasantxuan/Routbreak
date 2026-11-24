named_nb_matrix <- function(x, ids = paste0("c", seq_len(nrow(x)))) {
  dimnames(x) <- list(ids, ids)
  x
}

nb_reference_probability <- function(R, omega, tol = 1e-14,
                                     max_iter = 100000L) {
  q <- rep(0, nrow(R))
  for (iteration in seq_len(max_iter)) {
    q_next <- numeric(nrow(R))
    for (source in seq_len(ncol(R))) {
      log_q <- 0
      for (destination in seq_len(nrow(R))) {
        log_q <- log_q - log1p(
          omega * R[destination, source] * (1 - q[destination])
        ) / omega
      }
      q_next[source] <- exp(log_q)
    }
    if (max(abs(q_next - q)) <= tol) return(1 - q_next)
    q <- q_next
  }
  stop("Reference solver did not converge.")
}

nb_reference_rob <- function(R, p, omega) {
  ans <- numeric(ncol(R))
  for (source in seq_len(ncol(R))) {
    total <- 0
    for (destination in seq_len(nrow(R))) {
      total <- total +
        R[destination, source] * p[destination] *
        (1 + omega * R[destination, source]) /
        (1 + omega * R[destination, source] * p[destination])
    }
    ans[source] <- total / p[source]
  }
  ans
}

simulate_nb_survival <- function(R, source, omega, n_runs = 4000L,
                                 cutoff = 300L, max_generations = 100L) {
  survived <- logical(n_runs)
  for (run in seq_len(n_runs)) {
    infected <- integer(nrow(R))
    infected[source] <- 1L
    total_size <- 1L

    for (generation in seq_len(max_generations)) {
      next_infected <- integer(nrow(R))
      for (source_j in seq_len(ncol(R))) {
        if (infected[source_j] == 0L) next
        for (destination_i in seq_len(nrow(R))) {
          if (R[destination_i, source_j] == 0) next
          next_infected[destination_i] <-
            next_infected[destination_i] + rnbinom(
              1L,
              size = infected[source_j] / omega,
              mu = infected[source_j] * R[destination_i, source_j]
            )
        }
      }
      infected <- next_infected
      total_size <- total_size + sum(infected)
      if (sum(infected) == 0L) break
      if (total_size >= cutoff) {
        survived[run] <- TRUE
        break
      }
    }
  }
  mean(survived)
}

test_that("single-population negative-binomial benchmark is exact", {
  R <- matrix(2, 1, 1, dimnames = list("one", "one"))
  p <- outbreak_prob(R, offspring = "nb", omega = 1)
  expect_equal(p, c(one = 0.5), tolerance = 2e-12)
  expect_equal(r_ob(R, p, offspring = "nb", omega = 1),
               c(one = 3), tolerance = 1e-12)
})

test_that("negative-binomial calculations have the Poisson limit", {
  R <- named_nb_matrix(matrix(c(1.3, 0.05, 0.5, 0.8), 2, byrow = TRUE))
  poisson_p <- outbreak_prob(R)
  poisson_rob <- r_ob(R, poisson_p)

  expect_identical(
    outbreak_prob(R, offspring = "nb", omega = 0),
    poisson_p
  )
  expect_identical(
    r_ob(R, poisson_p, offspring = "nb", omega = 0),
    poisson_rob
  )
  expect_equal(
    outbreak_prob(R, offspring = "nb", omega = 1e-8),
    poisson_p,
    tolerance = 2e-7
  )
  expect_equal(
    r_ob(R, offspring = "nb", omega = 1e-8),
    poisson_rob,
    tolerance = 2e-7
  )
})

test_that("asymmetric negative-binomial orientation matches a reference solver", {
  R <- named_nb_matrix(matrix(c(1.3, 0.05, 0.5, 0.8), 2, byrow = TRUE))
  omega <- 0.7
  expected_p <- nb_reference_probability(R, omega)
  expected_rob <- nb_reference_rob(R, expected_p, omega)

  p <- outbreak_prob(R, offspring = "nb", omega = omega)
  rob <- r_ob(R, p, offspring = "nb", omega = omega)

  expect_equal(unname(p), expected_p, tolerance = 1e-10)
  expect_equal(unname(rob), expected_rob, tolerance = 1e-10)
  expect_equal(p, c(c1 = 0.3465643056036252, c2 = 0.0686317302829711),
               tolerance = 1e-10)
  expect_equal(rob, c(c1 = 2.018214452043478, c2 = 1.459996127389071),
               tolerance = 1e-10)
})

test_that("negative-binomial R_ob matches conditional expectations", {
  R <- matrix(c(1.3, 0.05, 0.5, 0.8), 2, byrow = TRUE)
  omega <- 0.7
  p <- outbreak_prob(R, offspring = "nb", omega = omega)
  q <- 1 - p
  counts <- 0:150

  enumerated <- vapply(seq_len(ncol(R)), function(source) {
    probability_1 <- dnbinom(
      counts, size = 1 / omega, mu = R[1, source]
    )
    probability_2 <- dnbinom(
      counts, size = 1 / omega, mu = R[2, source]
    )
    grid <- expand.grid(x1 = counts, x2 = counts)
    joint_probability <- as.vector(outer(probability_1, probability_2))
    total_offspring <- grid$x1 + grid$x2
    descendant_extinction <- q[1]^grid$x1 * q[2]^grid$x2

    epidemic_mean <- sum(
      total_offspring * joint_probability * (1 - descendant_extinction)
    ) / p[source]
    extinction_mean <- sum(
      total_offspring * joint_probability * descendant_extinction
    ) / q[source]
    epidemic_mean - extinction_mean
  }, numeric(1))

  expect_equal(
    unname(r_ob(R, p, offspring = "nb", omega = omega)),
    enumerated,
    tolerance = 1e-11
  )
})

test_that("negative-binomial probabilities are finite and bounded", {
  R <- named_nb_matrix(matrix(c(1.1, 0.2, 0.3, 1.0), 2, byrow = TRUE))
  for (omega in c(0.1, 1, 5)) {
    p <- outbreak_prob(R, offspring = "nb", omega = omega)
    expect_true(all(is.finite(p)))
    expect_true(all(p >= 0 & p <= 1))
    residual <- max(abs(
      p - vapply(seq_len(ncol(R)), function(source) {
        1 - prod((1 + omega * R[, source] * p)^(-1 / omega))
      }, numeric(1))
    ))
    expect_lt(residual, 1e-11)
  }
})

test_that("near-critical negative-binomial probabilities are stable", {
  subcritical <- matrix(0.99999, 1, 1,
                        dimnames = list("one", "one"))
  expect_equal(
    outbreak_prob(subcritical, max_iter = 1,
                  offspring = "nb", omega = 1),
    c(one = 0),
    tolerance = 0
  )

  mean <- 1.00001
  supercritical <- matrix(mean, 1, 1,
                          dimnames = list("one", "one"))
  expect_lt(
    abs(unname(outbreak_prob(
      supercritical, offspring = "nb", omega = 1
    )) - (mean - 1) / mean),
    1e-12
  )
})

test_that("negative-binomial validation and nonconvergence are explicit", {
  R <- matrix(2, 1, 1)
  expect_error(outbreak_prob(R, offspring = "geometric"), "offspring")
  expect_error(outbreak_prob(R, offspring = "nb", omega = -0.1),
               "nonnegative")
  expect_error(outbreak_prob(R, omega = 1), "must be zero")
  expect_error(
    outbreak_prob(R, max_iter = 1, offspring = "nb", omega = 1),
    "did not converge"
  )

  subcritical <- named_nb_matrix(
    matrix(c(0.4, 0.2, 0.05, 0.6), 2, byrow = TRUE)
  )
  expect_equal(
    outbreak_prob(subcritical, offspring = "nb", omega = 1),
    c(c1 = 0, c2 = 0),
    tolerance = 0
  )
  expect_error(
    r_ob(subcritical, offspring = "nb", omega = 1),
    "do not provide a negative-binomial continuation"
  )
})

test_that("seeded branching simulations validate NB outbreak probabilities", {
  R <- matrix(c(1.3, 0.05, 0.5, 0.8), 2, byrow = TRUE)
  omega <- 0.7
  theoretical <- outbreak_prob(R, offspring = "nb", omega = omega)

  set.seed(20261122)
  simulated <- vapply(seq_len(ncol(R)), function(source) {
    simulate_nb_survival(R, source, omega)
  }, numeric(1))

  expect_lt(max(abs(simulated - unname(theoretical))), 0.025)
})
