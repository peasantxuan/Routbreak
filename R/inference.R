#' @importFrom stats quantile rgamma setNames
NULL

# Validation helpers ------------------------------------------------------

.validate_count_value <- function(x, name, allow_zero = TRUE) {
  valid <- is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
    all(x == floor(x))
  if (allow_zero) valid <- valid && all(x >= 0) else valid <- valid && all(x > 0)
  if (!valid) {
    qualifier <- if (allow_zero) "nonnegative" else "positive"
    stop(sprintf("`%s` must contain finite %s integer counts.", name, qualifier),
         call. = FALSE)
  }
  invisible(x)
}

.validate_named_counts <- function(x, name) {
  .validate_count_value(x, name)
  if (is.null(names(x)) || anyNA(names(x)) || any(names(x) == "") ||
      anyDuplicated(names(x))) {
    stop(sprintf(
      "`%s` must have non-missing, non-empty, unique community identifiers.",
      name
    ), call. = FALSE)
  }
  invisible(x)
}

.expand_prior <- function(x, n, ids, name) {
  if (is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
      x > 0) {
    return(matrix(x, n, n, dimnames = list(ids, ids)))
  }
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(n, n)) ||
      any(!is.finite(x)) || any(x <= 0)) {
    stop(sprintf(
      "`%s` must be a finite positive scalar or a finite positive %d x %d matrix.",
      name, n, n
    ), call. = FALSE)
  }
  if (!identical(dimnames(x), list(ids, ids))) {
    stop(sprintf(
      "Matrix `%s` must have row and column identifiers matching `S` in order.",
      name
    ), call. = FALSE)
  }
  x
}

.validate_structural_zero <- function(x, n, ids) {
  if (is.null(x)) {
    return(matrix(FALSE, n, n, dimnames = list(ids, ids)))
  }
  if (!is.matrix(x) || !is.logical(x) || !identical(dim(x), c(n, n)) ||
      anyNA(x)) {
    stop(sprintf(
      "`structural_zero` must be NULL or a non-missing logical %d x %d matrix.",
      n, n
    ), call. = FALSE)
  }
  if (!identical(dimnames(x), list(ids, ids))) {
    stop(paste0(
      "`structural_zero` must have row and column identifiers matching `S` ",
      "in order."
    ), call. = FALSE)
  }
  x
}

.validate_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(x)
}

.validate_positive_integer <- function(x, name) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x >= 1 && x <= .Machine$integer.max && x == floor(x)
  if (!valid) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  as.integer(x)
}

.validate_seed <- function(x) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x >= 0 && x <= .Machine$integer.max && x == floor(x)
  if (!valid) {
    stop("`seed` must be a nonnegative integer.", call. = FALSE)
  }
  as.integer(x)
}

#' Infer a reproduction matrix from sufficient statistics
#'
#' Fits the independent Poisson transmission model in Eqs. (27)--(30) of
#' Wang and Valdano (2026). `S[i, j]` is the total number of observed
#' transmissions from source `j` to destination `i`, and `E[j]` is the total
#' number of infectious source-`j` individuals exposed across the preceding
#' generations. Rows are destinations and columns are sources.
#'
#' The maximum-likelihood estimate is `S[i, j] / E[j]` when `E[j] > 0`.
#' Entries with zero exposure are reported as unidentified rather than set to
#' zero. Under a proper Gamma prior in shape-rate form, posterior shape is
#' `prior_shape + S` and posterior rate is `prior_rate + E[j]`.
#'
#' Prior parameters are required explicitly because a fixed independent prior
#' on every edge implies a different prior for total reproduction as the
#' number of communities changes. For example, if every edge has a common
#' `Gamma(a, b)` prior, the prior mean of each source-column total is
#' `n * a / b`. Users should elicit priors for the matrix scale and structure
#' relevant to their study. This is not the paper's improper flat-prior
#' construction. The posterior mean is `(prior_shape + S) /
#' (prior_rate + E)`, whereas the MLE remains `S / E`.
#'
#' @param S A finite nonnegative square numeric matrix of integer transmission
#'   counts. Row and column identifiers must be unique and identical in the
#'   same order.
#' @param E A named finite nonnegative numeric vector of integer source
#'   exposures, in the same order as the columns of `S`.
#' @param prior_shape A required finite positive scalar or matrix of Gamma
#'   prior shape parameters.
#' @param prior_rate A required finite positive scalar or matrix of Gamma
#'   prior rate parameters.
#' @param structural_zero Optional logical matrix identifying edges fixed at
#'   zero by design. Positive observations on those edges are rejected. A zero
#'   observed count is otherwise treated as a sampling zero, not a structural
#'   zero.
#'
#' @return An object of class `routbreak_inference` containing sufficient
#'   statistics `S` and `E`; square destination-by-source matrices `mle`,
#'   `status`, `prior_shape`, `prior_rate`, `alpha_post`, `beta_post`,
#'   `posterior_mean`, and `posterior_variance`; the logical matrix
#'   `structural_zero`; and model metadata. The `mle` is `NA` for unidentified
#'   zero-exposure entries, while their proper-prior posterior moments remain
#'   defined. Declared structural zeros have zero MLE and posterior moments.
#' @examples
#' ids <- c("A", "B")
#' S <- matrix(c(12, 4, 3, 16), 2, byrow = TRUE,
#'             dimnames = list(ids, ids))
#' fit <- infer_R_sufficient(
#'   S, E = c(A = 10, B = 20), prior_shape = 1, prior_rate = 1
#' )
#' fit$mle
#' fit$posterior_mean
#' @export
infer_R_sufficient <- function(S, E, prior_shape = NULL,
                               prior_rate = NULL,
                               structural_zero = NULL) {
  .validate_R(S, require_ids = TRUE)
  .validate_count_value(as.numeric(S), "S")
  ids <- rownames(S)
  n <- nrow(S)

  if (is.null(prior_shape) || is.null(prior_rate)) {
    stop(paste0(
      "Specify both `prior_shape` and `prior_rate` explicitly. Routbreak ",
      "does not use a universal per-edge Gamma prior because its implied ",
      "total-reproduction prior depends on the matrix dimension."
    ), call. = FALSE)
  }

  .validate_named_counts(E, "E")
  if (!identical(names(E), ids) || length(E) != n) {
    stop("Names of `E` must match the matrix community identifiers in order.",
         call. = FALSE)
  }

  structural_zero <- .validate_structural_zero(structural_zero, n, ids)
  if (any(S[structural_zero] > 0)) {
    stop("Positive transmissions are incompatible with `structural_zero` edges.",
         call. = FALSE)
  }
  zero_exposure <- E == 0
  if (any(zero_exposure) && any(S[, zero_exposure, drop = FALSE] > 0)) {
    stop("Positive transmissions cannot arise from a source with zero exposure.",
         call. = FALSE)
  }

  prior_shape_matrix <- .expand_prior(prior_shape, n, ids, "prior_shape")
  prior_rate_matrix <- .expand_prior(prior_rate, n, ids, "prior_rate")
  exposure_matrix <- matrix(
    rep(as.numeric(E), each = n),
    nrow = n,
    dimnames = list(ids, ids)
  )

  mle <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  identified_columns <- E > 0
  if (any(identified_columns)) {
    mle[, identified_columns] <- sweep(
      S[, identified_columns, drop = FALSE],
      2L,
      E[identified_columns],
      `/`
    )
  }
  mle[structural_zero] <- 0

  status <- matrix(
    rep(ifelse(E > 0, "identified", "unidentified"), each = n),
    nrow = n,
    dimnames = list(ids, ids)
  )
  status[structural_zero] <- "structural_zero"

  alpha_post <- prior_shape_matrix + S
  beta_post <- prior_rate_matrix + exposure_matrix
  posterior_mean <- alpha_post / beta_post
  posterior_variance <- alpha_post / beta_post^2

  alpha_post[structural_zero] <- NA_real_
  beta_post[structural_zero] <- NA_real_
  posterior_mean[structural_zero] <- 0
  posterior_variance[structural_zero] <- 0

  structure(list(
    call = match.call(),
    model = "independent Poisson transmission with Gamma priors",
    data_mode = "sufficient-statistics",
    community = ids,
    S = S,
    E = E,
    mle = mle,
    status = status,
    structural_zero = structural_zero,
    prior_shape = prior_shape_matrix,
    prior_rate = prior_rate_matrix,
    alpha_post = alpha_post,
    beta_post = beta_post,
    posterior_mean = posterior_mean,
    posterior_variance = posterior_variance
  ), class = "routbreak_inference")
}

#' Infer a reproduction matrix from generation-resolved transmission data
#'
#' Converts generation-resolved directed transmission counts into the
#' sufficient statistics used by [infer_R_sufficient()]. For generation `t`,
#' `Y[t, i, j]` records infections in destination `i` caused by source `j`.
#' Source exposure at generation `t` is the number infected in source `j` at
#' generation `t - 1`; the initially infected generation is therefore included
#' in exposure for the first observed transmission generation.
#'
#' Missing source-destination rows are interpreted as zero transmissions.
#' Duplicate generation-source-destination rows are rejected rather than
#' silently aggregated.
#'
#' @param data A data frame with columns `generation`, `source`, `destination`,
#'   and `transmissions`. Generations are positive integers beginning with the
#'   first observed offspring generation.
#' @param initial_cases A named vector of finite nonnegative integer counts for
#'   generation zero. Its names define the canonical community order.
#' @inheritParams infer_R_sufficient
#'
#' @return An object of class `routbreak_inference`. In addition to the fields
#'   returned by [infer_R_sufficient()], it contains the reconstructed
#'   array `Y` with dimensions generation by destination by source, named
#'   `initial_cases`, and the generation-by-community matrix
#'   `infections_by_generation` used to calculate previous-generation
#'   exposures.
#' @examples
#' transmissions <- data.frame(
#'   generation = rep(1L, 4),
#'   source = c("A", "B", "A", "B"),
#'   destination = c("A", "A", "B", "B"),
#'   transmissions = c(12L, 4L, 3L, 16L)
#' )
#' fit <- infer_R(
#'   transmissions, initial_cases = c(A = 10, B = 20),
#'   prior_shape = 1, prior_rate = 1
#' )
#' fit$S
#' fit$E
#' fit$mle
#' @export
infer_R <- function(data, initial_cases, prior_shape = NULL,
                    prior_rate = NULL, structural_zero = NULL) {
  required <- c("generation", "source", "destination", "transmissions")
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop(sprintf(
      "`data` must be a data frame containing columns: %s.",
      paste(required, collapse = ", ")
    ), call. = FALSE)
  }
  if (!nrow(data)) stop("`data` must contain at least one observation.", call. = FALSE)

  .validate_named_counts(initial_cases, "initial_cases")
  ids <- names(initial_cases)
  n <- length(ids)

  .validate_count_value(data$generation, "data$generation", allow_zero = FALSE)
  .validate_count_value(data$transmissions, "data$transmissions")
  if (!is.character(data$source) || !is.character(data$destination) ||
      anyNA(data$source) || anyNA(data$destination) ||
      any(data$source == "") || any(data$destination == "")) {
    stop("Transmission source and destination identifiers must be non-missing character values.",
         call. = FALSE)
  }
  unknown <- setdiff(unique(c(data$source, data$destination)), ids)
  if (length(unknown)) {
    stop(sprintf(
      "Transmission data contain communities absent from `initial_cases`: %s.",
      paste(unknown, collapse = ", ")
    ), call. = FALSE)
  }

  keys <- paste(data$generation, data$source, data$destination, sep = "\r")
  if (anyDuplicated(keys)) {
    stop("Transmission data contain duplicate generation-source-destination rows.",
         call. = FALSE)
  }

  generations <- as.integer(data$generation)
  T_max <- max(generations)
  Y <- array(
    0,
    dim = c(T_max, n, n),
    dimnames = list(
      generation = as.character(seq_len(T_max)),
      destination = ids,
      source = ids
    )
  )
  for (row in seq_len(nrow(data))) {
    Y[generations[row], data$destination[row], data$source[row]] <-
      data$transmissions[row]
  }

  infections <- matrix(
    0,
    nrow = T_max + 1L,
    ncol = n,
    dimnames = list(generation = as.character(0:T_max), community = ids)
  )
  infections[1L, ] <- as.numeric(initial_cases)

  for (generation in seq_len(T_max)) {
    exposure <- infections[generation, ]
    transmissions_by_source <- apply(
      Y[generation, , , drop = FALSE], 3L, sum
    )
    if (any(exposure == 0 & transmissions_by_source > 0)) {
      bad <- ids[exposure == 0 & transmissions_by_source > 0]
      stop(sprintf(
        paste0("Generation %d contains positive transmissions from source(s) ",
               "with zero previous-generation exposure: %s."),
        generation, paste(bad, collapse = ", ")
      ), call. = FALSE)
    }
    infections[generation + 1L, ] <- apply(
      Y[generation, , , drop = FALSE], 2L, sum
    )
  }

  S <- apply(Y, c(2L, 3L), sum)
  dimnames(S) <- list(ids, ids)
  E <- setNames(colSums(infections[seq_len(T_max), , drop = FALSE]), ids)

  fit <- infer_R_sufficient(
    S = S,
    E = E,
    prior_shape = prior_shape,
    prior_rate = prior_rate,
    structural_zero = structural_zero
  )
  fit$call <- match.call()
  fit$data_mode <- "generation-resolved"
  fit$initial_cases <- initial_cases
  fit$Y <- Y
  fit$infections_by_generation <- infections
  fit
}

#' Propagate reproduction-matrix uncertainty to outbreak risk
#'
#' Draws complete reproduction matrices from a fitted Gamma posterior and
#' applies [outbreak_prob()], [r_ob()], and optionally [r_ref()] to every draw.
#' The transmission likelihood is Poisson; negative-binomial transmission
#' inference is not implied by the offspring-distribution option in
#' [outbreak_prob()]. Intervals are posterior credible intervals, not
#' frequentist bootstrap intervals.
#'
#' Draws that fail numerical or theoretical checks are recorded. Any failures
#' produce a warning, and the function errors when their fraction exceeds
#' `max_fail_fraction`. R's Gamma generator can underflow to exact zero for
#' very small shape parameters even though a proper Gamma distribution has
#' strictly positive support. Such values are replaced by the smallest
#' positive normalized floating-point value and counted in the diagnostics;
#' this prevents numerical underflow from creating artificial structural
#' zeros.
#'
#' Posterior draws are independent Gamma samples, not an MCMC chain, so MCMC
#' convergence diagnostics do not apply. Reported quantiles still have Monte
#' Carlo error; users should increase `draws` and check stability for final
#' analyses. If some draws fail within the allowed fraction, summaries are
#' conditional on the successful draws and must be interpreted with the
#' accompanying failure diagnostics.
#'
#' @param fit An object returned by [infer_R()] or [infer_R_sufficient()].
#' @param draws Positive integer number of posterior matrices to sample.
#' @param seed Optional finite integer random seed.
#' @param level Credible interval level strictly between zero and one.
#' @param return_draws If `TRUE`, include posterior matrix, outbreak-ratio, and
#'   outbreak-probability draws in the returned object.
#' @param include_rref If `TRUE`, calculate the reference reproduction ratio
#'   for every successful draw.
#' @param max_fail_fraction Maximum tolerated fraction of failed draws.
#' @param tol,max_iter Numerical controls passed to [outbreak_prob()].
#'
#' @return An object of class `routbreak_rob_estimate`. Its `summary` data
#'   frame contains posterior medians and lower/upper credible limits for
#'   `R_ob`, outbreak probability, and (when requested) `R_ref`.
#'   `requested_draws`, `successful_draws`, `failed_draws`, `failures`, and
#'   `diagnostics` report the failure fraction, reasons, numerical controls,
#'   and Gamma-underflow replacements. If `return_draws = TRUE`, successful
#'   matrix, `R_ob`, outbreak-probability, and optional `R_ref` draws are
#'   returned in arrays or matrices preserving community identifiers.
#' @examples
#' S <- matrix(12, 1, 1, dimnames = list("A", "A"))
#' fit <- infer_R_sufficient(
#'   S, E = c(A = 10), prior_shape = 1, prior_rate = 1
#' )
#' estimate <- estimate_rob(fit, draws = 20, seed = 1)
#' estimate$summary
#' estimate$diagnostics
#' @export
estimate_rob <- function(fit, draws = 1000, seed = NULL, level = 0.95,
                         return_draws = FALSE, include_rref = TRUE,
                         max_fail_fraction = 0.05, tol = 1e-12,
                         max_iter = 100000) {
  if (!inherits(fit, "routbreak_inference")) {
    stop("`fit` must be an object returned by `infer_R()` or `infer_R_sufficient()`.",
         call. = FALSE)
  }
  draws <- .validate_positive_integer(draws, "draws")
  .validate_flag(return_draws, "return_draws")
  .validate_flag(include_rref, "include_rref")
  .validate_positive_scalar(tol, "tol")
  max_iter <- .validate_positive_integer(max_iter, "max_iter")

  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a finite numeric scalar strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.numeric(max_fail_fraction) || length(max_fail_fraction) != 1L ||
      !is.finite(max_fail_fraction) || max_fail_fraction < 0 ||
      max_fail_fraction >= 1) {
    stop("`max_fail_fraction` must be in [0, 1).", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .validate_seed(seed)
    had_random_seed <- exists(".Random.seed", envir = .GlobalEnv,
                              inherits = FALSE)
    if (had_random_seed) {
      caller_random_seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
    }
    on.exit({
      if (had_random_seed) {
        assign(".Random.seed", caller_random_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  ids <- fit$community
  n <- length(ids)
  matrix_draws <- array(
    0,
    dim = c(draws, n, n),
    dimnames = list(draw = as.character(seq_len(draws)),
                    destination = ids, source = ids)
  )
  rob_draws <- matrix(NA_real_, draws, n,
                      dimnames = list(draw = as.character(seq_len(draws)),
                                      community = ids))
  probability_draws <- rob_draws
  rref_draws <- rep(NA_real_, draws)
  failures <- vector("list", draws)
  underflow_replacements <- 0L

  sampled_edges <- !fit$structural_zero
  for (draw in seq_len(draws)) {
    sampled <- matrix(0, n, n, dimnames = list(ids, ids))
    sampled_values <- rgamma(
      sum(sampled_edges),
      shape = fit$alpha_post[sampled_edges],
      rate = fit$beta_post[sampled_edges]
    )
    underflow <- sampled_values == 0
    underflow_replacements <- underflow_replacements + sum(underflow)
    sampled_values[underflow] <- .Machine$double.xmin
    sampled[sampled_edges] <- sampled_values
    matrix_draws[draw, , ] <- sampled

    result <- tryCatch({
      p <- outbreak_prob(sampled, tol = tol, max_iter = max_iter)
      rob <- r_ob(sampled, p)
      list(p = p, rob = rob,
           rref = if (include_rref) r_ref(sampled) else NA_real_)
    }, error = function(e) e)

    if (inherits(result, "error")) {
      failures[[draw]] <- conditionMessage(result)
    } else {
      probability_draws[draw, ] <- result$p
      rob_draws[draw, ] <- result$rob
      rref_draws[draw] <- result$rref
    }
  }

  failed_index <- which(lengths(failures) > 0L)
  successful_index <- setdiff(seq_len(draws), failed_index)
  failure_fraction <- length(failed_index) / draws
  if (failure_fraction > max_fail_fraction) {
    reasons <- sort(table(unlist(failures[failed_index])), decreasing = TRUE)
    stop(sprintf(
      paste0("%d of %d posterior draws failed (%.1f%%), exceeding ",
             "`max_fail_fraction = %.3g`. Most common reason: %s"),
      length(failed_index), draws, 100 * failure_fraction,
      max_fail_fraction, names(reasons)[1L]
    ), call. = FALSE)
  }
  if (!length(successful_index)) {
    stop("No posterior draw produced a supported outbreak-risk calculation.",
         call. = FALSE)
  }
  if (length(failed_index)) {
    warning(sprintf(
      "%d of %d posterior draws failed; inspect `$failures` for reasons.",
      length(failed_index), draws
    ), call. = FALSE)
  }

  alpha <- (1 - level) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  summarize_columns <- function(x) {
    t(apply(x[successful_index, , drop = FALSE], 2L, quantile,
            probs = probs, names = FALSE))
  }
  rob_summary <- summarize_columns(rob_draws)
  p_summary <- summarize_columns(probability_draws)
  if (include_rref) {
    rref_summary <- quantile(rref_draws[successful_index], probs = probs,
                             names = FALSE)
  } else {
    rref_summary <- rep(NA_real_, 3L)
  }

  summary <- data.frame(
    community = ids,
    R_ob_median = rob_summary[, 2L],
    R_ob_lower = rob_summary[, 1L],
    R_ob_upper = rob_summary[, 3L],
    outbreak_probability_median = p_summary[, 2L],
    outbreak_probability_lower = p_summary[, 1L],
    outbreak_probability_upper = p_summary[, 3L],
    R_ref_median = rep(rref_summary[2L], n),
    R_ref_lower = rep(rref_summary[1L], n),
    R_ref_upper = rep(rref_summary[3L], n),
    stringsAsFactors = FALSE
  )

  failure_table <- if (length(failed_index)) {
    data.frame(
      draw = failed_index,
      reason = unlist(failures[failed_index]),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(draw = integer(), reason = character(),
               stringsAsFactors = FALSE)
  }

  structure(list(
    call = match.call(),
    summary = summary,
    level = level,
    requested_draws = draws,
    successful_draws = length(successful_index),
    failed_draws = length(failed_index),
    failures = failure_table,
    diagnostics = list(
      tolerance = tol,
      max_iterations = max_iter,
      failure_fraction = failure_fraction,
      gamma_underflow_replacements = underflow_replacements,
      convergence = paste0(
        "Every successful draw satisfied the fixed-point residual check ",
        "implemented by outbreak_prob()."
      )
    ),
    seed = seed,
    include_rref = include_rref,
    successful_draw_indices = successful_index,
    matrix_draws = if (return_draws) {
      matrix_draws[successful_index, , , drop = FALSE]
    } else NULL,
    R_ob_draws = if (return_draws) rob_draws[successful_index, , drop = FALSE] else NULL,
    outbreak_probability_draws = if (return_draws) {
      probability_draws[successful_index, , drop = FALSE]
    } else NULL,
    R_ref_draws = if (return_draws && include_rref) {
      rref_draws[successful_index]
    } else NULL
  ), class = "routbreak_rob_estimate")
}
