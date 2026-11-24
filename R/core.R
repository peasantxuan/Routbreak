# Internal validation -----------------------------------------------------

.validate_R <- function(R, require_ids = FALSE) {
  if (!is.matrix(R) || !is.numeric(R)) {
    stop("`R` must be a numeric matrix.", call. = FALSE)
  }
  if (length(dim(R)) != 2L || nrow(R) < 1L || nrow(R) != ncol(R)) {
    stop("`R` must be a non-empty square matrix.", call. = FALSE)
  }
  if (any(!is.finite(R))) {
    stop("`R` must contain only finite values; missing values are not allowed.",
         call. = FALSE)
  }
  if (any(R < 0)) {
    stop("`R` must be nonnegative.", call. = FALSE)
  }

  rn <- rownames(R)
  cn <- colnames(R)
  if (is.null(rn) != is.null(cn)) {
    stop("`R` must have both row and column identifiers, or neither.",
         call. = FALSE)
  }
  if (!is.null(rn)) {
    invalid <- anyNA(rn) || anyNA(cn) || any(rn == "") || any(cn == "") ||
      anyDuplicated(rn) || anyDuplicated(cn)
    if (invalid) {
      stop("Matrix community identifiers must be non-missing, non-empty, and unique.",
           call. = FALSE)
    }
    if (!identical(rn, cn)) {
      stop(paste0(
        "Row and column community identifiers must match in the same order; ",
        "Routbreak does not silently reorder matrices."
      ), call. = FALSE)
    }
  } else if (require_ids) {
    stop("The matrix must have matching row and column community identifiers.",
         call. = FALSE)
  }
  invisible(R)
}

.community_names <- function(R) {
  if (is.null(colnames(R))) as.character(seq_len(ncol(R))) else colnames(R)
}

.spectral_radius <- function(R) {
  max(Mod(eigen(R, only.values = TRUE)$values))
}

.is_diagonal <- function(R) {
  n <- nrow(R)
  if (n == 1L) return(TRUE)
  all(R[row(R) != col(R)] == 0)
}

.is_irreducible <- function(R) {
  n <- nrow(R)
  if (n == 1L) return(TRUE)
  reachable <- R > 0
  diag(reachable) <- TRUE
  for (k in seq_len(n)) {
    reachable <- reachable | (reachable[, k] %o% reachable[k, ])
  }
  all(reachable)
}

.strong_components <- function(R) {
  n <- nrow(R)
  reachable <- R > 0
  diag(reachable) <- TRUE
  for (k in seq_len(n)) {
    reachable <- reachable | (reachable[, k] %o% reachable[k, ])
  }

  component <- integer(n)
  next_component <- 0L
  for (seed in seq_len(n)) {
    if (component[seed] != 0L) next
    next_component <- next_component + 1L
    members <- which(reachable[, seed] & reachable[seed, ])
    component[members] <- next_component
  }
  component
}

.integer_gcd <- function(a, b) {
  a <- abs(as.integer(a))
  b <- abs(as.integer(b))
  while (b != 0L) {
    remainder <- a %% b
    a <- b
    b <- remainder
  }
  a
}

.component_period <- function(R, members) {
  adjacency <- t(R) > 0
  depth <- rep(NA_integer_, nrow(R))
  root <- members[1L]
  depth[root] <- 0L
  queue <- root
  head <- 1L
  while (head <= length(queue)) {
    source <- queue[head]
    head <- head + 1L
    destinations <- members[adjacency[source, members]]
    new_destinations <- destinations[is.na(depth[destinations])]
    if (length(new_destinations)) {
      depth[new_destinations] <- depth[source] + 1L
      queue <- c(queue, new_destinations)
    }
  }

  edges <- which(adjacency[members, members, drop = FALSE], arr.ind = TRUE)
  if (!nrow(edges)) return(NA_integer_)
  period <- 0L
  for (edge in seq_len(nrow(edges))) {
    source <- members[edges[edge, 1L]]
    destination <- members[edges[edge, 2L]]
    difference <- depth[source] + 1L - depth[destination]
    period <- .integer_gcd(period, difference)
  }
  if (period == 0L) NA_integer_ else period
}

.matrix_period <- function(R) {
  .component_period(R, seq_len(nrow(R)))
}

.validate_positive_scalar <- function(x, name, allow_zero = FALSE) {
  valid_scalar <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  lower_ok <- valid_scalar && if (allow_zero) x >= 0 else x > 0
  if (!lower_ok) {
    qualifier <- if (allow_zero) "nonnegative" else "positive"
    stop(sprintf("`%s` must be a finite %s numeric scalar.", name, qualifier),
         call. = FALSE)
  }
  invisible(x)
}

.validate_offspring <- function(offspring, omega) {
  if (!is.character(offspring) || length(offspring) != 1L ||
      is.na(offspring) || !offspring %in% c("poisson", "nb")) {
    stop("`offspring` must be either \"poisson\" or \"nb\".", call. = FALSE)
  }
  .validate_positive_scalar(omega, "omega", allow_zero = TRUE)
  if (offspring == "poisson" && omega != 0) {
    stop("`omega` must be zero when `offspring = \"poisson\"`.", call. = FALSE)
  }
  if (offspring == "nb" && omega == 0) offspring <- "poisson"
  list(offspring = offspring, omega = omega)
}

.validate_supplied_probability <- function(R, p, offspring, omega) {
  nm <- .community_names(R)
  if (!is.numeric(p) || is.matrix(p) || length(p) != nrow(R) ||
      any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop(paste0("`p` must be a finite numeric vector of length nrow(R), ",
                "with every value between zero and one."), call. = FALSE)
  }
  if (!is.null(names(p)) && !identical(names(p), nm)) {
    stop("Names of `p` must match the matrix community identifiers in order.",
         call. = FALSE)
  }

  expected <- outbreak_prob(R, offspring = offspring, omega = omega)
  validation_tol <- max(100 * .Machine$double.eps, 100 * 1e-12)
  if (max(abs(as.numeric(p) - expected)) > validation_tol) {
    stop(paste0(
      "Supplied `p` does not match the maximal outbreak-probability ",
      "fixed point for `R` within the validation tolerance."
    ), call. = FALSE)
  }
  ans <- as.numeric(p)
  names(ans) <- nm
  ans
}

#' Outbreak probabilities for a spatial branching process
#'
#' Computes the maximal fixed point of the multitype branching-process
#' equations in Wang and Valdano (2026), using the convention that rows of
#' `R` are destinations and columns are sources. Iteration starts at a vector
#' of ones. The Poisson calculation uses
#' `-expm1(-crossprod(R, p))`. The negative-binomial calculation implements
#' Eq. (5) with `log1p()` and `expm1()` for numerical stability.
#'
#' In subcritical and critical regimes, the theoretically correct maximal
#' fixed point is zero. The function returns it directly when the spectral
#' radius is at or below one (within a relative numerical tolerance), avoiding
#' arbitrarily slow convergence near the threshold. Above the threshold,
#' monotone iteration is supplemented by safeguarded Newton steps when useful;
#' the fixed-point residual is still checked before a result is returned.
#'
#' @param R A finite, nonnegative square reproduction matrix. `R[i, j]` is
#'   the expected number of infections in destination `i` caused by one
#'   infected individual from source `j`.
#' @param tol Positive convergence tolerance for both the iteration change
#'   and fixed-point residual.
#' @param max_iter Positive integer giving the maximum number of iterations.
#' @param offspring Offspring distribution. Either `"poisson"` (the default)
#'   or `"nb"` for independent negative-binomial counts on directed edges.
#' @param omega A finite nonnegative scalar overdispersion parameter. For the
#'   negative-binomial model, `size = 1 / omega` and
#'   `Var(X[i,j]) = R[i,j] * (1 + omega * R[i,j])`. The value zero recovers
#'   the Poisson model. Heterogeneous dispersion matrices are unsupported.
#'
#' @return A named numeric vector, indexed by source community.
#' @examples
#' R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
#'             dimnames = list(c("A", "B"), c("A", "B")))
#' outbreak_prob(R)
#' @export
outbreak_prob <- function(R, tol = 1e-12, max_iter = 100000,
                          offspring = "poisson", omega = 0) {
  .validate_R(R)
  .validate_positive_scalar(tol, "tol")
  model <- .validate_offspring(offspring, omega)
  valid_max_iter <- is.numeric(max_iter) && length(max_iter) == 1L &&
    !is.na(max_iter) && is.finite(max_iter) && max_iter >= 1 &&
    max_iter <= .Machine$integer.max && max_iter == floor(max_iter)
  if (!valid_max_iter) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  n <- nrow(R)
  nm <- .community_names(R)
  rho <- .spectral_radius(R)
  critical_tol <- 100 * .Machine$double.eps * max(1, rho)

  if (rho <= 1 + critical_tol) {
    ans <- rep(0, n)
    names(ans) <- nm
    return(ans)
  }

  if (model$offspring == "poisson") {
    update <- function(p) -expm1(-as.vector(crossprod(R, p)))
    jacobian <- function(p) {
      extinction <- exp(-as.vector(crossprod(R, p)))
      sweep(t(R), 1L, extinction, `*`)
    }
  } else {
    update <- function(p) {
      edge_risk <- sweep(R, 1L, p, `*`)
      log_extinction <- -colSums(log1p(model$omega * edge_risk)) /
        model$omega
      -expm1(log_extinction)
    }
    jacobian <- function(p) {
      edge_risk <- sweep(R, 1L, p, `*`)
      log_extinction <- -colSums(log1p(model$omega * edge_risk)) /
        model$omega
      edge_derivative <- R / (1 + model$omega * edge_risk)
      sweep(t(edge_derivative), 1L, exp(log_extinction), `*`)
    }
  }
  p <- rep(1, n)
  converged <- FALSE
  residual <- Inf
  estimated_error <- Inf

  for (iteration in seq_len(max_iter)) {
    p_next <- update(p)
    change <- max(abs(p_next - p))
    residual <- max(abs(update(p_next) - p_next))
    p <- p_next
    if (change <= tol && residual <= tol) {
      local_correction <- tryCatch(
        solve(diag(n) - jacobian(p), update(p) - p),
        error = function(e) NULL
      )
      estimated_error <- if (is.null(local_correction) ||
                             any(!is.finite(local_correction))) {
        change
      } else {
        max(abs(local_correction))
      }
      if (estimated_error <= tol) {
        converged <- TRUE
        break
      }
    }

    # Fixed-point convergence becomes arbitrarily slow immediately above the
    # threshold. A safeguarded Newton step accelerates the same fixed-point
    # equation without accepting values outside the probability cube or a
    # step that increases the residual. The following ordinary iteration
    # still confirms convergence before return.
    if (iteration >= 100L && iteration < max_iter) {
      fixed_point_error <- update(p) - p
      newton_step <- tryCatch(
        solve(diag(n) - jacobian(p), fixed_point_error),
        error = function(e) NULL
      )
      if (!is.null(newton_step) && all(is.finite(newton_step))) {
        step_size <- 1
        for (line_search in seq_len(20L)) {
          candidate <- p + step_size * newton_step
          valid_candidate <- all(is.finite(candidate)) &&
            all(candidate >= 0 & candidate <= 1) && any(candidate > 0)
          if (valid_candidate) {
            candidate_residual <- max(abs(update(candidate) - candidate))
            if (candidate_residual < residual) {
              p <- candidate
              residual <- candidate_residual
              break
            }
          }
          step_size <- step_size / 2
        }
      }
    }
  }

  if (!converged) {
    stop(sprintf(
      paste0("Outbreak-probability iteration did not converge after %d ",
             "iterations (fixed-point residual %.3g; estimated solution ",
             "error %.3g; tolerance %.3g)."),
      max_iter, residual, estimated_error, tol
    ), call. = FALSE)
  }

  names(p) <- nm
  p
}

#' Outbreak reproduction ratio
#'
#' For Poisson offspring, implements Eq. (3) and the below-threshold
#' continuation in Eqs. (20)--(26) of Wang and Valdano (2026). For
#' supercritical irreducible matrices, the stable expression
#' `crossprod(R, p) / p` is used. For subcritical or critical irreducible
#' primitive matrices, the Perron--Frobenius continuation is the spectral
#' radius for every community. An irreducible but periodic matrix can lack the
#' required generation-by-generation limit and is reported explicitly.
#' Diagonal matrices reduce community-by-community to `diag(R)`.
#'
#' For independent negative-binomial counts on directed edges, implements
#' Eq. (S4) of the S1 Appendix. This expression is derived directly from the
#' conditional-expectation definition of the outbreak reproduction ratio.
#' The publication does not establish a negative-binomial continuation when
#' an outbreak probability is zero, so that regime is rejected explicitly.
#'
#' General reducible, non-diagonal matrices are unsupported for the Poisson
#' continuation because the published result assumes irreducibility; no
#' component convention is imposed by this package.
#'
#' @inheritParams outbreak_prob
#' @param p Optional outbreak-probability vector. If supplied, it must match
#'   the maximal fixed-point solution and, when named, the matrix identifiers.
#'
#' @return A named numeric vector of outbreak reproduction ratios.
#' @examples
#' R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
#'             dimnames = list(c("A", "B"), c("A", "B")))
#' r_ob(R)
#' @export
r_ob <- function(R, p = NULL, offspring = "poisson", omega = 0) {
  .validate_R(R)
  model <- .validate_offspring(offspring, omega)

  if (is.null(p)) {
    p <- outbreak_prob(R, offspring = offspring, omega = omega)
  } else {
    p <- .validate_supplied_probability(R, p, offspring, omega)
  }

  .r_ob_from_validated_p(R, p, model)
}

.r_ob_from_validated_p <- function(R, p, model) {
  nm <- .community_names(R)

  if (model$offspring == "nb") {
    if (any(p <= 0)) {
      stop(paste0(
        "`r_ob()` for negative-binomial offspring requires positive outbreak ",
        "probabilities. Wang and Valdano (2026) do not provide a ",
        "negative-binomial continuation for zero-probability source ",
        "communities."
      ), call. = FALSE)
    }
    edge_risk <- sweep(R, 1L, p, `*`)
    terms <- edge_risk * (1 + model$omega * R) /
      (1 + model$omega * edge_risk)
    ans <- colSums(terms) / p
  } else if (.is_diagonal(R)) {
    ans <- diag(R)
  } else {
    if (!.is_irreducible(R)) {
      stop(paste0(
        "`r_ob()` does not support reducible non-diagonal matrices: the ",
        "Perron--Frobenius continuation in Wang and Valdano (2026) assumes ",
        "a strongly connected (irreducible) reproduction matrix."
      ), call. = FALSE)
    }
    rho <- .spectral_radius(R)
    critical_tol <- 100 * .Machine$double.eps * max(1, rho)
    if (rho <= 1 + critical_tol) {
      period <- .matrix_period(R)
      if (!is.na(period) && period > 1L) {
        stop(paste0(
          "The Poisson continuation is not defined by a single generation-",
          "limit for this irreducible but periodic matrix (cyclic period ",
          period, "). Inspect `matrix_diagnostics(R)`; Routbreak does not ",
          "replace the published limit with an arbitrary cycle average."
        ), call. = FALSE)
      }
      ans <- rep(rho, nrow(R))
    } else {
      if (any(p <= 0)) {
        stop("A supercritical irreducible matrix must have positive outbreak probabilities.",
             call. = FALSE)
      }
      ans <- as.vector(crossprod(R, p)) / p
    }
  }
  names(ans) <- nm
  ans
}

#' Summarize outbreak risk by community
#'
#' Calculates outbreak probabilities once and reuses the validated result to
#' report the principal community-level quantities. Rows of `R` are
#' destinations and columns are sources.
#'
#' @inheritParams outbreak_prob
#' @return A data frame with one row per source community and columns
#'   `community`, `outbreak_probability`, `R_ob`, `R_local`, `R_total`, and
#'   `R_ref`.
#' @examples
#' R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
#'             dimnames = list(c("A", "B"), c("A", "B")))
#' rob_table(R)
#' @export
rob_table <- function(R, offspring = "poisson", omega = 0) {
  .validate_R(R)
  model <- .validate_offspring(offspring, omega)
  p <- outbreak_prob(R, offspring = offspring, omega = omega)
  rob <- .r_ob_from_validated_p(R, p, model)
  data.frame(
    community = .community_names(R),
    outbreak_probability = as.numeric(p),
    R_ob = as.numeric(rob),
    R_local = as.numeric(r_local(R)),
    R_total = as.numeric(r_total(R)),
    R_ref = rep(r_ref(R), nrow(R)),
    stringsAsFactors = FALSE
  )
}

#' Convert directed long-format data to a reproduction matrix
#'
#' Creates a square matrix whose rows are destinations and columns are
#' sources. An explicit zero-valued row creates a zero matrix entry. An absent
#' source-destination pair remains an unobserved edge and causes an error by
#' default; it is converted to zero only when `missing = "zero"` is requested.
#'
#' @param data A data frame containing directed edge data.
#' @param source,destination,value Single strings naming the source,
#'   destination, and nonnegative numeric value columns.
#' @param communities Optional character vector defining the complete set and
#'   order of community identifiers. When omitted, order is determined by the
#'   first appearance among destinations and then sources.
#' @param missing How absent directed edges are handled. The default,
#'   `"error"`, rejects incomplete edge data. Use `"zero"` only when absence
#'   is scientifically intended to mean a zero matrix entry.
#'
#' @return A numeric reproduction matrix with destination rows, source
#'   columns, and named `destination` and `source` dimname axes.
#' @examples
#' edges <- data.frame(
#'   source = c("A", "A", "B", "B"),
#'   destination = c("A", "B", "A", "B"),
#'   value = c(1.2, 0.3, 0.2, 0.8)
#' )
#' as_rmatrix(edges, "source", "destination", "value")
#' @export
as_rmatrix <- function(data, source, destination, value,
                       communities = NULL, missing = c("error", "zero")) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  column_args <- list(source = source, destination = destination, value = value)
  valid_names <- vapply(column_args, function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }, logical(1))
  if (!all(valid_names)) {
    stop("`source`, `destination`, and `value` must be single column-name strings.",
         call. = FALSE)
  }
  requested <- unlist(column_args, use.names = FALSE)
  if (anyDuplicated(requested)) {
    stop("`source`, `destination`, and `value` must name different columns.",
         call. = FALSE)
  }
  absent_columns <- setdiff(requested, names(data))
  if (length(absent_columns)) {
    stop(sprintf("Missing data column(s): %s.",
                 paste(absent_columns, collapse = ", ")), call. = FALSE)
  }
  if (!nrow(data)) stop("`data` must contain at least one directed edge.", call. = FALSE)
  missing <- match.arg(missing)

  source_id <- data[[source]]
  destination_id <- data[[destination]]
  if (is.factor(source_id)) source_id <- as.character(source_id)
  if (is.factor(destination_id)) destination_id <- as.character(destination_id)
  if (!is.character(source_id) || !is.character(destination_id) ||
      anyNA(source_id) || anyNA(destination_id) ||
      any(source_id == "") || any(destination_id == "")) {
    stop("Source and destination identifiers must be non-missing, non-empty character values.",
         call. = FALSE)
  }
  weights <- data[[value]]
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("Edge values must be finite, nonnegative numeric values.", call. = FALSE)
  }

  if (is.null(communities)) {
    communities <- unique(c(destination_id, source_id))
  } else if (!is.character(communities) || !length(communities) ||
             anyNA(communities) || any(communities == "") ||
             anyDuplicated(communities)) {
    stop("`communities` must contain non-missing, non-empty, unique character identifiers.",
         call. = FALSE)
  }
  unknown <- setdiff(unique(c(source_id, destination_id)), communities)
  if (length(unknown)) {
    stop(sprintf("Edge data contain unknown communities: %s.",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }

  keys <- paste(destination_id, source_id, sep = "\r")
  if (anyDuplicated(keys)) {
    stop("Edge data contain duplicate destination-source pairs.", call. = FALSE)
  }
  expected_count <- length(communities)^2
  if (length(keys) != expected_count && missing == "error") {
    stop(sprintf(
      paste0("Directed edge data are incomplete: found %d of %d destination-",
             "source pairs. Use `missing = \"zero\"` only when absent edges ",
             "are scientifically intended to be zero."),
      length(keys), expected_count
    ), call. = FALSE)
  }

  ans <- matrix(
    0, length(communities), length(communities),
    dimnames = list(destination = communities, source = communities)
  )
  ans[cbind(match(destination_id, communities),
            match(source_id, communities))] <- weights
  ans
}

#' Test whether a reproduction matrix is irreducible
#'
#' A matrix is irreducible when its directed transmission graph is strongly
#' connected: every community can reach every other community through one or
#' more positive-weight edges.
#'
#' @inheritParams outbreak_prob
#' @return A single logical value.
#' @examples
#' R <- matrix(c(1, 0.2, 0.1, 1), 2, byrow = TRUE,
#'             dimnames = list(c("A", "B"), c("A", "B")))
#' is_irreducible(R)
#' @export
is_irreducible <- function(R) {
  .validate_R(R)
  .is_irreducible(R)
}

#' Diagnose connected components and periodicity
#'
#' Reports strongly connected components of the directed reproduction graph.
#' `cyclic_period = 1` denotes a primitive component; values greater than one
#' denote periodic (imprimitive) components. For an isolated vertex without a
#' self-loop, the cyclic period is undefined and reported as `NA`.
#'
#' Periodicity is separate from reducibility. In particular, the published
#' generation-limit continuation used by [r_ob()] need not exist for an
#' irreducible but periodic subcritical or critical matrix. Routbreak reports
#' that case rather than replacing the limit with a cycle average.
#'
#' @inheritParams outbreak_prob
#' @return A data frame with one row per community and its component number,
#'   component size, cyclic period, primitive status, and matrix-level
#'   irreducibility status.
#' @examples
#' R <- matrix(c(1, 0.2, 0.1, 1), 2, byrow = TRUE,
#'             dimnames = list(c("A", "B"), c("A", "B")))
#' matrix_diagnostics(R)
#' @export
matrix_diagnostics <- function(R) {
  .validate_R(R)
  ids <- .community_names(R)
  component <- .strong_components(R)
  component_ids <- sort(unique(component))
  sizes <- tabulate(component, nbins = length(component_ids))
  periods <- vapply(component_ids, function(id) {
    .component_period(R, which(component == id))
  }, integer(1))
  data.frame(
    community = ids,
    component = component,
    component_size = sizes[component],
    cyclic_period = periods[component],
    primitive_component = !is.na(periods[component]) & periods[component] == 1L,
    matrix_irreducible = rep(length(component_ids) == 1L, length(ids)),
    stringsAsFactors = FALSE
  )
}

#' Reference reproduction ratio
#'
#' @inheritParams outbreak_prob
#' @return A single numeric value equal to the spectral radius of `R`.
#' @examples
#' R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE)
#' r_ref(R)
#' @export
r_ref <- function(R) {
  .validate_R(R)
  .spectral_radius(R)
}

#' Local reproduction ratios
#'
#' @inheritParams outbreak_prob
#' @return A named numeric vector equal to `diag(R)`.
#' @examples
#' R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
#'             dimnames = list(c("A", "B"), c("A", "B")))
#' r_local(R)
#' @export
r_local <- function(R) {
  .validate_R(R)
  ans <- diag(R)
  names(ans) <- .community_names(R)
  ans
}

#' Total reproduction ratios
#'
#' Returns the expected total number of infections generated anywhere by an
#' infected individual from each source community.
#'
#' @inheritParams outbreak_prob
#' @return A named numeric vector equal to `colSums(R)`.
#' @examples
#' R <- matrix(c(1.2, 0.2, 0.3, 0.8), 2, byrow = TRUE,
#'             dimnames = list(c("A", "B"), c("A", "B")))
#' r_total(R)
#' @export
r_total <- function(R) {
  .validate_R(R)
  ans <- colSums(R)
  names(ans) <- .community_names(R)
  ans
}

#' Scale a reproduction matrix to a target reference ratio
#'
#' @inheritParams outbreak_prob
#' @param target A finite nonnegative target spectral radius.
#' @return A scaled reproduction matrix with the original dimnames.
#' @examples
#' R <- matrix(c(1, 0.2, 0.1, 0.8), 2, byrow = TRUE)
#' scaled <- scale_to_rref(R, target = 1.5)
#' r_ref(scaled)
#' @export
scale_to_rref <- function(R, target) {
  .validate_R(R)
  .validate_positive_scalar(target, "target", allow_zero = TRUE)
  rho <- .spectral_radius(R)
  if (rho == 0) {
    if (target > 0) {
      stop("Cannot scale a zero-spectral-radius matrix to a positive target.",
           call. = FALSE)
    }
    return(R)
  }
  R * (target / rho)
}

#' Construct a reproduction matrix from colocation data
#'
#' Constructs `A[i, j] = C[i, j] * population[i]`, then scales `A` so its
#' spectral radius is `target_rref`. Community identifiers must be explicit,
#' unique, and identical across matrix rows, matrix columns, and the population
#' vector; inputs are never silently reordered.
#'
#' @param C A finite, nonnegative square numeric colocation matrix with
#'   matching row and column identifiers.
#' @param population A named finite positive numeric vector of destination
#'   population sizes.
#' @param target_rref A finite positive target reference reproduction ratio.
#'
#' @return A reproduction matrix with the same dimnames as `C`.
#' @examples
#' C <- matrix(c(0.8, 0.1, 0.2, 0.7), 2, byrow = TRUE,
#'             dimnames = list(c("A", "B"), c("A", "B")))
#' population <- c(A = 1000, B = 2000)
#' R <- colocation_to_r(C, population, target_rref = 1.5)
#' r_ref(R)
#' @export
colocation_to_r <- function(C, population, target_rref) {
  .validate_R(C, require_ids = TRUE)
  .validate_positive_scalar(target_rref, "target_rref")
  ids <- rownames(C)
  if (!is.numeric(population) || is.matrix(population) ||
      length(population) != nrow(C) || any(!is.finite(population)) ||
      any(population <= 0)) {
    stop(paste0("`population` must be a finite positive numeric vector with ",
                "one value per community."), call. = FALSE)
  }
  if (is.null(names(population)) || anyNA(names(population)) ||
      any(names(population) == "") || anyDuplicated(names(population))) {
    stop("`population` must have non-missing, non-empty, unique community identifiers.",
         call. = FALSE)
  }
  if (!identical(names(population), ids)) {
    stop(paste0(
      "Population identifiers must match matrix identifiers in the same ",
      "order; Routbreak does not silently reorder data."
    ), call. = FALSE)
  }

  A <- C * population
  dimnames(A) <- dimnames(C)
  rho <- .spectral_radius(A)
  if (rho == 0) {
    stop("The population-weighted colocation matrix has zero spectral radius and cannot be scaled to a positive target.",
         call. = FALSE)
  }
  A * (target_rref / rho)
}
