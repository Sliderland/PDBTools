# Posterior prediction and diagnostics for Bayesian VECM(p) models.
#
# Model convention:
#   Delta y[t] = intercept + Pi %*% y[t - 1] +
#                sum_{i=1}^{p-1} Gamma_i %*% Delta y[t - i] + error[t]
#   Pi = alpha %*% t(beta)
#
# Array conventions:
# - alpha_draws: draws x variables x rank
# - beta_draws:  draws x variables x rank
# - gamma_draws: draws x (p - 1) x variables x variables. For p = 1,
#                supply an array with zero in its second dimension.
# - sigma_draws: draws x variables x variables
# - fitted/forecast draws: draws x time-or-horizon x variables
#
# Typical use with an rstan fit:
#   posterior <- extract_vecm_draws(fit)
#   fitted <- vecm_fitted_draws(
#       posterior$alpha_draws, posterior$beta_draws,
#       posterior$gamma_draws, y
#   )
#   forecasts <- vecm_forecast_draws(
#       posterior$alpha_draws, posterior$beta_draws,
#       posterior$gamma_draws, posterior$sigma_draws,
#       y, horizon = 8L
#   )
#   residual_checks <- posterior_vecm_residual_diagnostics(
#       posterior$alpha_draws, posterior$beta_draws,
#       posterior$gamma_draws, y
#   )
#   rank_and_roots <- vecm_rank_root_diagnostics(
#       posterior$alpha_draws, posterior$beta_draws,
#       posterior$gamma_draws
#   )

validate_vecm_inputs <- function(alpha_draws, beta_draws, gamma_draws, y,
                                 sigma_draws = NULL, intercept = NULL) {
    y <- as.matrix(y)
    alpha_draws <- as.array(alpha_draws)
    beta_draws <- as.array(beta_draws)
    gamma_draws <- as.array(gamma_draws)

    if (length(dim(alpha_draws)) != 3L ||
        length(dim(beta_draws)) != 3L) {
        stop(
            "`alpha_draws` and `beta_draws` must have dimensions draws x variables x rank.",
            call. = FALSE
        )
    }
    if (length(dim(gamma_draws)) != 4L) {
        stop(
            "`gamma_draws` must have dimensions draws x (p - 1) x variables x variables.",
            call. = FALSE
        )
    }
    if (!is.numeric(y) || !is.numeric(alpha_draws) ||
        !is.numeric(beta_draws) || !is.numeric(gamma_draws)) {
        stop("All data and parameter arrays must be numeric.", call. = FALSE)
    }

    n_draws <- dim(alpha_draws)[1L]
    k <- ncol(y)
    rank <- dim(alpha_draws)[3L]
    q <- dim(gamma_draws)[2L]
    p <- q + 1L

    if (rank < 1L || rank > k) {
        stop("The cointegration rank must lie between 1 and ncol(y).", call. = FALSE)
    }
    if (!identical(dim(alpha_draws), c(n_draws, k, rank)) ||
        !identical(dim(beta_draws), c(n_draws, k, rank))) {
        stop("`alpha_draws` and `beta_draws` have incompatible dimensions.", call. = FALSE)
    }
    if (!identical(dim(gamma_draws), c(n_draws, q, k, k))) {
        stop("`gamma_draws` has incompatible dimensions.", call. = FALSE)
    }
    if (nrow(y) <= p) {
        stop("`y` must contain more rows than the level-VAR lag order.", call. = FALSE)
    }

    finite_objects <- list(y, alpha_draws, beta_draws, gamma_draws)
    if (any(vapply(finite_objects, function(x) {
        anyNA(x) || any(!is.finite(x))
    }, logical(1)))) {
        stop("Data and parameter arrays must contain only finite values.", call. = FALSE)
    }

    if (!is.null(sigma_draws)) {
        sigma_draws <- as.array(sigma_draws)
        if (!identical(dim(sigma_draws), c(n_draws, k, k))) {
            stop(
                "`sigma_draws` must have dimensions draws x variables x variables.",
                call. = FALSE
            )
        }
        if (anyNA(sigma_draws) || any(!is.finite(sigma_draws))) {
            stop("`sigma_draws` must contain only finite values.", call. = FALSE)
        }
        for (d in seq_len(n_draws)) {
            if (!isTRUE(all.equal(
                sigma_draws[d, , ], t(sigma_draws[d, , ]), tolerance = 1e-10
            ))) {
                stop("Every covariance draw must be symmetric.", call. = FALSE)
            }
            if (inherits(try(chol(sigma_draws[d, , ]), silent = TRUE), "try-error")) {
                stop("Every covariance draw must be positive definite.", call. = FALSE)
            }
        }
    }

    if (is.null(intercept)) {
        intercept <- numeric(k)
    }
    if (is.atomic(intercept) && is.null(dim(intercept)) && length(intercept) == k) {
        intercept <- matrix(
            rep(intercept, times = n_draws), nrow = n_draws, byrow = TRUE
        )
    } else {
        intercept <- as.matrix(intercept)
    }
    if (!identical(dim(intercept), c(n_draws, k)) ||
        anyNA(intercept) || any(!is.finite(intercept))) {
        stop(
            "`intercept` must be a finite length-k vector or draws x k matrix.",
            call. = FALSE
        )
    }

    list(
        alpha_draws = alpha_draws,
        beta_draws = beta_draws,
        gamma_draws = gamma_draws,
        sigma_draws = sigma_draws,
        y = y,
        intercept = intercept,
        n_draws = n_draws,
        rank = rank,
        q = q,
        p = p,
        k = k
    )
}

select_vecm_draws <- function(n_draws, draws = NULL, max_draws = NULL) {
    if (is.null(draws)) {
        draws <- seq_len(n_draws)
    }
    if (!is.numeric(draws) || anyNA(draws) || any(draws %% 1 != 0) ||
        any(draws < 1L | draws > n_draws)) {
        stop("`draws` contains invalid draw indices.", call. = FALSE)
    }
    draws <- unique(as.integer(draws))
    if (!is.null(max_draws)) {
        if (length(max_draws) != 1L || !is.finite(max_draws) ||
            max_draws < 1L || max_draws %% 1 != 0) {
            stop("`max_draws` must be a positive integer.", call. = FALSE)
        }
        if (length(draws) > max_draws) {
            draws <- draws[unique(round(seq(1, length(draws), length.out = max_draws)))]
        }
    }
    draws
}

vecm_posterior_point_estimate <- function(draws, method = c("median", "mean")) {
    method <- match.arg(method)
    draws <- as.array(draws)
    if (length(dim(draws)) < 2L) {
        stop("`draws` must have a leading posterior-draw dimension.", call. = FALSE)
    }
    margin <- seq.int(2L, length(dim(draws)))
    if (method == "median") {
        apply(draws, margin, stats::median)
    } else {
        apply(draws, margin, mean)
    }
}

extract_vecm_draws <- function(x, alpha_name = "alpha", beta_name = "beta",
                               gamma_name = "Gamma", sigma_name = "Sigma",
                               intercept_name = NULL, include_sigma = TRUE) {
    parameters <- c(
        alpha_name, beta_name, gamma_name,
        if (include_sigma) sigma_name,
        if (!is.null(intercept_name)) intercept_name
    )
    extracted <- if (inherits(x, "stanfit")) {
        if (!requireNamespace("rstan", quietly = TRUE)) {
            stop("The `rstan` package is required for a `stanfit` input.", call. = FALSE)
        }
        rstan::extract(x, pars = parameters, permuted = TRUE)
    } else if (is.list(x)) {
        x
    } else {
        stop("`x` must be a `stanfit` or an extracted-draw list.", call. = FALSE)
    }

    required <- c(alpha_name, beta_name, gamma_name, if (include_sigma) sigma_name)
    missing <- required[vapply(required, function(name) {
        is.null(extracted[[name]])
    }, logical(1))]
    if (length(missing)) {
        stop("No draws found for: ", paste(missing, collapse = ", "), ".", call. = FALSE)
    }

    result <- list(
        alpha_draws = extracted[[alpha_name]],
        beta_draws = extracted[[beta_name]],
        gamma_draws = extracted[[gamma_name]]
    )
    if (include_sigma) {
        result$sigma_draws <- extracted[[sigma_name]]
    }
    if (!is.null(intercept_name)) {
        if (is.null(extracted[[intercept_name]])) {
            stop("No intercept draws named `", intercept_name, "` were found.", call. = FALSE)
        }
        result$intercept_draws <- extracted[[intercept_name]]
    }
    result
}

# Extract the parameterization used by the VECM_info.RDS archive created by
# unconstrainedVECM.R. The saved Stan models use beta as rank x variables,
# xi as the short-run matrices, L as an innovation Cholesky factor, and mu as
# the unrestricted intercept. This adapter returns the conventions expected by
# the remaining functions in this file.
extract_saved_vecm_draws <- function(archive, model_name,
                                     use_saved_draws = TRUE) {
    if (is.character(archive) && length(archive) == 1L) {
        archive <- readRDS(archive)
    }
    if (!is.list(archive) || is.null(archive$models) || is.null(archive$samps)) {
        stop(
            "`archive` must be a VECM_info.RDS path or its deserialized list.",
            call. = FALSE
        )
    }
    available <- intersect(names(archive$models), names(archive$samps))
    if (length(model_name) != 1L || !model_name %in% available) {
        stop(
            "Unknown `model_name`. Available models: ",
            paste(available, collapse = ", "), ".",
            call. = FALSE
        )
    }
    if (!requireNamespace("posterior", quietly = TRUE)) {
        stop("The `posterior` package is required for saved CmdStan draws.", call. = FALSE)
    }

    draws <- if (use_saved_draws) {
        archive$samps[[model_name]]
    } else {
        archive$models[[model_name]]$draws()
    }
    variables <- posterior::variables(draws)
    has_variable <- function(name) {
        any(variables == name | startsWith(variables, paste0(name, "[")))
    }
    short_run_name <- if (has_variable("xi")) {
        "xi"
    } else if (has_variable("gamma")) {
        "gamma"
    } else {
        "xi"
    }
    required <- c(
        "alpha", "beta", short_run_name, "L", "mu",
        if (has_variable("phi")) "phi"
    )
    missing <- required[!vapply(required, has_variable, logical(1))]
    if (length(missing)) {
        stop(
            "The saved draws do not contain: ", paste(missing, collapse = ", "),
            ".", call. = FALSE
        )
    }
    # Subset first: the archive also contains large transformed/generated
    # quantities, and converting all of them would need substantial memory.
    draws <- posterior::subset_draws(draws, variable = required)
    rvars <- posterior::as_draws_rvars(draws)

    flatten_rvar <- function(x) {
        values <- posterior::draws_of(x)
        dimensions <- dim(values)
        # as_draws_rvars() has already combined iteration and chain into the
        # leading draw dimension.
        array(values, dim = c(dimensions[1L], dimensions[-1L]))
    }

    alpha <- flatten_rvar(rvars$alpha)
    beta_stan <- flatten_rvar(rvars$beta)
    short_run <- flatten_rvar(rvars[[short_run_name]])
    chol_draws <- flatten_rvar(rvars$L)
    mu <- flatten_rvar(rvars$mu)

    n_draws <- dim(alpha)[1L]
    k <- dim(alpha)[2L]
    rank <- dim(alpha)[3L]
    if (!identical(dim(beta_stan), c(n_draws, rank, k))) {
        stop("The saved `beta` dimensions do not match `alpha`.", call. = FALSE)
    }
    beta <- array(NA_real_, dim = c(n_draws, k, rank))
    sigma <- array(NA_real_, dim = c(n_draws, k, k))
    for (d in seq_len(n_draws)) {
        beta[d, , ] <- t(matrix(beta_stan[d, , ], nrow = rank, ncol = k))
        chol <- matrix(chol_draws[d, , ], nrow = k, ncol = k)
        sigma[d, , ] <- tcrossprod(chol)
    }

    result <- list(
        alpha_draws = alpha,
        beta_draws = beta,
        gamma_draws = short_run,
        sigma_draws = sigma,
        intercept_draws = matrix(mu, nrow = n_draws, ncol = k),
        model_name = model_name
    )
    if (has_variable("phi")) {
        result$exogenous_draws <- flatten_rvar(rvars$phi)
    }
    if (model_name %in% c("vecm_long_run", "vecm_urca_hmc")) {
        extra_note <- if (identical(model_name, "vecm_urca_hmc")) {
            " It also includes exogenous centered seasonal indicators."
        } else {
            ""
        }
        warning(
            "`", model_name,
            "` uses y[t - p] in its error-correction term; ",
            "the fitted, residual, forecast, and VECM-to-VAR functions in ",
            "this file use the transitory y[t - 1] convention.",
            extra_note,
            call. = FALSE
        )
    }
    result
}

vecm_pi_draws <- function(alpha_draws, beta_draws) {
    alpha_draws <- as.array(alpha_draws)
    beta_draws <- as.array(beta_draws)
    if (length(dim(alpha_draws)) != 3L ||
        !identical(dim(alpha_draws), dim(beta_draws))) {
        stop("`alpha_draws` and `beta_draws` must have identical 3D dimensions.", call. = FALSE)
    }
    out <- array(NA_real_, dim = c(dim(alpha_draws)[1L], dim(alpha_draws)[2L],
                                   dim(alpha_draws)[2L]))
    for (d in seq_len(dim(out)[1L])) {
        alpha <- matrix(alpha_draws[d, , ], nrow = dim(alpha_draws)[2L])
        beta <- matrix(beta_draws[d, , ], nrow = dim(beta_draws)[2L])
        out[d, , ] <- alpha %*% t(beta)
    }
    out
}

vecm_fitted_draws <- function(alpha_draws, beta_draws, gamma_draws, y,
                              intercept = NULL, draws = NULL,
                              max_draws = NULL) {
    inputs <- validate_vecm_inputs(
        alpha_draws, beta_draws, gamma_draws, y, intercept = intercept
    )
    selected <- select_vecm_draws(inputs$n_draws, draws, max_draws)
    dy <- diff(inputs$y)
    fitted <- array(
        NA_real_,
        dim = c(length(selected), nrow(inputs$y) - inputs$p, inputs$k),
        dimnames = list(
            draw = selected,
            time = seq.int(inputs$p + 1L, nrow(inputs$y)),
            variable = colnames(inputs$y)
        )
    )

    for (d_out in seq_along(selected)) {
        d <- selected[[d_out]]
        alpha <- matrix(inputs$alpha_draws[d, , ], nrow = inputs$k)
        beta <- matrix(inputs$beta_draws[d, , ], nrow = inputs$k)
        pi <- alpha %*% t(beta)
        for (t in seq.int(inputs$p + 1L, nrow(inputs$y))) {
            conditional_mean <- inputs$intercept[d, ] + pi %*% inputs$y[t - 1L, ]
            if (inputs$q > 0L) {
                for (lag in seq_len(inputs$q)) {
                    conditional_mean <- conditional_mean +
                        inputs$gamma_draws[d, lag, , ] %*% dy[t - lag - 1L, ]
                }
            }
            fitted[d_out, t - inputs$p, ] <- conditional_mean
        }
    }
    fitted
}

vecm_fitted_point <- function(alpha_draws, beta_draws, gamma_draws, y,
                              intercept = NULL,
                              method = c("median", "mean")) {
    method <- match.arg(method)
    alpha <- vecm_posterior_point_estimate(alpha_draws, method)
    beta <- vecm_posterior_point_estimate(beta_draws, method)
    gamma <- vecm_posterior_point_estimate(gamma_draws, method)
    alpha <- array(alpha, dim = c(1L, dim(alpha)))
    beta <- array(beta, dim = c(1L, dim(beta)))
    gamma <- array(gamma, dim = c(1L, dim(gamma)))
    point_intercept <- if (is.null(intercept) ||
        is.atomic(intercept) && is.null(dim(intercept))) {
        intercept
    } else {
        vecm_posterior_point_estimate(intercept, method)
    }
    drop(vecm_fitted_draws(
        alpha, beta, gamma, y, intercept = point_intercept
    )[1L, , , drop = FALSE])
}

vecm_residual_draws <- function(alpha_draws, beta_draws, gamma_draws, y,
                                intercept = NULL, draws = NULL,
                                max_draws = NULL) {
    fitted <- vecm_fitted_draws(
        alpha_draws, beta_draws, gamma_draws, y,
        intercept = intercept, draws = draws, max_draws = max_draws
    )
    p <- dim(as.array(gamma_draws))[2L] + 1L
    observed <- diff(as.matrix(y))[-seq_len(p - 1L), , drop = FALSE]
    sweep(fitted, c(2L, 3L), observed, "-") * -1
}

vecm_rmvn_chol <- function(sigma) {
    drop(t(chol(sigma)) %*% stats::rnorm(nrow(sigma)))
}

vecm_forecast_draws <- function(alpha_draws, beta_draws, gamma_draws,
                                sigma_draws, y, horizon, intercept = NULL,
                                draws = NULL, max_draws = NULL, seed = NULL) {
    inputs <- validate_vecm_inputs(
        alpha_draws, beta_draws, gamma_draws, y,
        sigma_draws = sigma_draws, intercept = intercept
    )
    if (length(horizon) != 1L || !is.finite(horizon) || horizon < 1L ||
        horizon %% 1 != 0) {
        stop("`horizon` must be a positive integer.", call. = FALSE)
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    selected <- select_vecm_draws(inputs$n_draws, draws, max_draws)
    forecasts <- array(
        NA_real_, dim = c(length(selected), as.integer(horizon), inputs$k),
        dimnames = list(
            draw = selected, horizon = seq_len(horizon), variable = colnames(inputs$y)
        )
    )

    for (d_out in seq_along(selected)) {
        d <- selected[[d_out]]
        history <- inputs$y
        alpha <- matrix(inputs$alpha_draws[d, , ], nrow = inputs$k)
        beta <- matrix(inputs$beta_draws[d, , ], nrow = inputs$k)
        pi <- alpha %*% t(beta)
        for (step in seq_len(horizon)) {
            dy_history <- diff(history)
            conditional_mean <- inputs$intercept[d, ] +
                pi %*% history[nrow(history), ]
            if (inputs$q > 0L) {
                for (lag in seq_len(inputs$q)) {
                    conditional_mean <- conditional_mean +
                        inputs$gamma_draws[d, lag, , ] %*%
                            dy_history[nrow(dy_history) - lag + 1L, ]
                }
            }
            next_difference <- drop(conditional_mean) +
                vecm_rmvn_chol(inputs$sigma_draws[d, , ])
            next_level <- history[nrow(history), ] + next_difference
            forecasts[d_out, step, ] <- next_level
            history <- rbind(history, next_level)
        }
    }
    forecasts
}

summarize_vecm_posterior_array <- function(x, probs = c(0.05, 0.5, 0.95)) {
    x <- as.array(x)
    if (length(dim(x)) != 3L) {
        stop("`x` must have dimensions draws x index x variables.", call. = FALSE)
    }
    if (any(probs < 0 | probs > 1)) {
        stop("Every value in `probs` must lie in [0, 1].", call. = FALSE)
    }
    list(
        mean = apply(x, c(2L, 3L), mean),
        quantiles = apply(x, c(2L, 3L), stats::quantile, probs = probs),
        probs = probs
    )
}

vecm_residual_acf <- function(residuals, lag_max = 20L,
                              correlation = TRUE) {
    residuals <- as.matrix(residuals)
    if (nrow(residuals) <= lag_max) {
        stop("`lag_max` must be smaller than the number of residuals.", call. = FALSE)
    }
    marginal <- lapply(seq_len(ncol(residuals)), function(j) {
        drop(stats::acf(
            residuals[, j], lag.max = lag_max, plot = FALSE,
            na.action = stats::na.fail
        )$acf)
    })
    names(marginal) <- colnames(residuals)
    cross <- lapply(seq_len(lag_max), function(lag) {
        stats::cov(
            residuals[(lag + 1L):nrow(residuals), , drop = FALSE],
            residuals[seq_len(nrow(residuals) - lag), , drop = FALSE]
        )
    })
    if (correlation) {
        scales <- sqrt(diag(stats::cov(residuals)))
        cross <- lapply(cross, function(x) x / outer(scales, scales))
    }
    list(marginal = marginal, cross = cross, lag_max = lag_max)
}

vecm_multivariate_portmanteau <- function(residuals, lags = 16L,
                                          fitted_difference_lags = 0L,
                                          adjusted = TRUE) {
    input_name <- deparse(substitute(residuals))
    residuals <- as.matrix(residuals)
    n <- nrow(residuals)
    k <- ncol(residuals)
    if (lags <= fitted_difference_lags || lags >= n) {
        stop(
            "Require `fitted_difference_lags < lags < nrow(residuals)`.",
            call. = FALSE
        )
    }
    c0 <- crossprod(residuals) / n
    c0_inverse <- tryCatch(
        solve(c0),
        error = function(error) stop("Residual covariance is singular.", call. = FALSE)
    )
    terms <- vapply(seq_len(lags), function(lag) {
        cj <- crossprod(
            residuals[(lag + 1L):n, , drop = FALSE],
            residuals[seq_len(n - lag), , drop = FALSE]
        ) / n
        sum(diag(t(cj) %*% c0_inverse %*% cj %*% c0_inverse))
    }, numeric(1))
    statistic <- if (adjusted) {
        n^2 * sum(terms / (n - seq_len(lags)))
    } else {
        n * sum(terms)
    }
    degrees_freedom <- k^2 * (lags - fitted_difference_lags)
    structure(
        list(
            statistic = c("Chi-squared" = statistic),
            parameter = c(df = degrees_freedom),
            p.value = stats::pchisq(statistic, degrees_freedom, lower.tail = FALSE),
            method = paste(
                if (adjusted) "Adjusted" else "Asymptotic",
                "multivariate Portmanteau test for VECM residuals"
            ),
            data.name = input_name
        ),
        class = "htest"
    )
}

vecm_vech_rows <- function(residuals) {
    residuals <- as.matrix(residuals)
    keep <- lower.tri(matrix(0, ncol(residuals), ncol(residuals)), diag = TRUE)
    t(vapply(seq_len(nrow(residuals)), function(i) {
        tcrossprod(residuals[i, ])[keep]
    }, numeric(sum(keep))))
}

vecm_multivariate_arch_lm <- function(residuals, lags = 1L) {
    input_name <- deparse(substitute(residuals))
    residuals <- scale(as.matrix(residuals))
    n <- nrow(residuals)
    k <- ncol(residuals)
    if (lags < 1L || lags >= n) {
        stop("Require `1 <= lags < nrow(residuals)`.", call. = FALSE)
    }
    squared_products <- vecm_vech_rows(residuals)
    embedded <- embed(squared_products, lags + 1L)
    n_products <- ncol(squared_products)
    response <- embedded[, seq_len(n_products), drop = FALSE]
    predictors <- embedded[, -seq_len(n_products), drop = FALSE]
    if (nrow(response) <= ncol(predictors) + 1L) {
        stop("The multivariate ARCH regression is underdetermined.", call. = FALSE)
    }
    unrestricted <- stats::lm.fit(cbind(1, predictors), response)$residuals
    restricted <- stats::lm.fit(matrix(1, nrow(response), 1L), response)$residuals
    omega1 <- stats::cov(unrestricted)
    omega0 <- stats::cov(restricted)
    r_squared <- 1 - (2 / (k * (k + 1))) *
        sum(diag(omega1 %*% solve(omega0)))
    statistic <- 0.5 * nrow(unrestricted) * k * (k + 1) * r_squared
    degrees_freedom <- lags * k^2 * (k + 1)^2 / 4
    structure(
        list(
            statistic = c("Chi-squared" = statistic),
            parameter = c(df = degrees_freedom),
            p.value = stats::pchisq(statistic, degrees_freedom, lower.tail = FALSE),
            r.squared = r_squared,
            method = "Multivariate ARCH-LM test for VECM residuals",
            data.name = input_name
        ),
        class = "htest"
    )
}

vecm_cointegration_projection_draws <- function(beta_draws, tolerance = 1e-10) {
    beta_draws <- as.array(beta_draws)
    if (length(dim(beta_draws)) != 3L) {
        stop("`beta_draws` must have dimensions draws x variables x rank.", call. = FALSE)
    }
    n_draws <- dim(beta_draws)[1L]
    k <- dim(beta_draws)[2L]
    out <- array(NA_real_, dim = c(n_draws, k, k))
    for (d in seq_len(n_draws)) {
        beta <- matrix(beta_draws[d, , ], nrow = k)
        gram <- crossprod(beta)
        if (min(svd(gram, nu = 0L, nv = 0L)$d) <= tolerance) {
            stop("A `beta` draw is not full column rank.", call. = FALSE)
        }
        out[d, , ] <- beta %*% solve(gram, t(beta))
    }
    out
}

vecm_to_var_draws <- function(alpha_draws, beta_draws, gamma_draws) {
    alpha_draws <- as.array(alpha_draws)
    beta_draws <- as.array(beta_draws)
    gamma_draws <- as.array(gamma_draws)
    if (length(dim(alpha_draws)) != 3L ||
        !identical(dim(alpha_draws), dim(beta_draws)) ||
        length(dim(gamma_draws)) != 4L ||
        dim(gamma_draws)[1L] != dim(alpha_draws)[1L] ||
        any(dim(gamma_draws)[3:4] != dim(alpha_draws)[2L])) {
        stop("The VECM parameter arrays have incompatible dimensions.", call. = FALSE)
    }
    n_draws <- dim(alpha_draws)[1L]
    k <- dim(alpha_draws)[2L]
    q <- dim(gamma_draws)[2L]
    p <- q + 1L
    phi <- array(0, dim = c(n_draws, p, k, k))
    identity <- diag(k)
    for (d in seq_len(n_draws)) {
        alpha <- matrix(alpha_draws[d, , ], nrow = k)
        beta <- matrix(beta_draws[d, , ], nrow = k)
        pi <- alpha %*% t(beta)
        if (q == 0L) {
            phi[d, 1L, , ] <- identity + pi
        } else {
            phi[d, 1L, , ] <- identity + pi + gamma_draws[d, 1L, , ]
            if (q > 1L) {
                for (lag in 2:q) {
                    phi[d, lag, , ] <- gamma_draws[d, lag, , ] -
                        gamma_draws[d, lag - 1L, , ]
                }
            }
            phi[d, p, , ] <- -gamma_draws[d, q, , ]
        }
    }
    phi
}

vecm_companion_roots <- function(alpha_draws, beta_draws, gamma_draws) {
    phi <- vecm_to_var_draws(alpha_draws, beta_draws, gamma_draws)
    n_draws <- dim(phi)[1L]
    p <- dim(phi)[2L]
    k <- dim(phi)[3L]
    roots <- matrix(NA_complex_, nrow = n_draws, ncol = k * p)
    for (d in seq_len(n_draws)) {
        companion <- matrix(0, nrow = k * p, ncol = k * p)
        companion[seq_len(k), ] <- do.call(cbind, lapply(seq_len(p), function(lag) {
            phi[d, lag, , ]
        }))
        if (p > 1L) {
            companion[(k + 1L):(k * p), seq_len(k * (p - 1L))] <-
                diag(k * (p - 1L))
        }
        roots[d, ] <- eigen(companion, only.values = TRUE)$values
    }
    roots
}

vecm_rank_root_diagnostics <- function(alpha_draws, beta_draws, gamma_draws,
                                       rank_tolerance = 1e-8,
                                       unit_root_tolerance = 1e-6) {
    pi_draws <- vecm_pi_draws(alpha_draws, beta_draws)
    roots <- vecm_companion_roots(alpha_draws, beta_draws, gamma_draws)
    expected_rank <- dim(as.array(alpha_draws))[3L]
    k <- dim(as.array(alpha_draws))[2L]
    numerical_rank <- vapply(seq_len(dim(pi_draws)[1L]), function(d) {
        singular_values <- svd(pi_draws[d, , ], nu = 0L, nv = 0L)$d
        sum(singular_values > rank_tolerance * max(1, singular_values[[1L]]))
    }, integer(1))
    unit_roots <- rowSums(abs(Mod(roots) - 1) <= unit_root_tolerance)
    explosive_roots <- rowSums(Mod(roots) > 1 + unit_root_tolerance)
    data.frame(
        draw = seq_len(nrow(roots)),
        numerical_rank = numerical_rank,
        expected_rank = expected_rank,
        unit_roots = unit_roots,
        expected_unit_roots = k - expected_rank,
        explosive_roots = explosive_roots,
        rank_ok = numerical_rank == expected_rank,
        unit_roots_ok = unit_roots == k - expected_rank,
        no_explosive_roots = explosive_roots == 0L
    )
}

posterior_vecm_residual_diagnostics <- function(
    alpha_draws,
    beta_draws,
    gamma_draws,
    y,
    intercept = NULL,
    draws = NULL,
    max_draws = 200L,
    serial_lags = 16L,
    arch_lags = 1L,
    alpha_level = 0.05
) {
    residual_draws <- vecm_residual_draws(
        alpha_draws, beta_draws, gamma_draws, y,
        intercept = intercept, draws = draws, max_draws = max_draws
    )
    q <- dim(as.array(gamma_draws))[2L]
    diagnostics <- lapply(seq_len(dim(residual_draws)[1L]), function(d) {
        list(
            serial = vecm_multivariate_portmanteau(
                residual_draws[d, , ], lags = serial_lags,
                fitted_difference_lags = q, adjusted = TRUE
            ),
            arch = vecm_multivariate_arch_lm(
                residual_draws[d, , ], lags = arch_lags
            )
        )
    })
    serial_p <- vapply(diagnostics, function(x) x$serial$p.value, numeric(1))
    arch_p <- vapply(diagnostics, function(x) x$arch$p.value, numeric(1))
    list(
        draw_indices = as.integer(dimnames(residual_draws)[[1L]]),
        residual_draws = residual_draws,
        diagnostics = diagnostics,
        summary = data.frame(
            test = c("serial_correlation", "multivariate_arch"),
            median_p_value = c(stats::median(serial_p), stats::median(arch_p)),
            proportion_rejecting = c(
                mean(serial_p < alpha_level), mean(arch_p < alpha_level)
            ),
            alpha = alpha_level,
            row.names = NULL
        )
    )
}

classical_vecm_diagnostics <- function(y, p, rank, deterministic = "const",
                                       serial_lags = 16L,
                                       arch_lags = 5L) {
    if (!requireNamespace("urca", quietly = TRUE) ||
        !requireNamespace("vars", quietly = TRUE)) {
        stop(
            "Install the `urca` and `vars` packages to use `classical_vecm_diagnostics()`.",
            call. = FALSE
        )
    }
    if (rank < 1L || rank >= ncol(as.matrix(y))) {
        stop("`rank` must lie between 1 and ncol(y) - 1.", call. = FALSE)
    }
    johansen <- urca::ca.jo(
        as.matrix(y), K = p, ecdet = deterministic, type = "trace"
    )
    vecm <- urca::cajorls(johansen, r = rank)
    level_var <- vars::vec2var(johansen, r = rank)
    list(
        johansen = johansen,
        vecm = vecm,
        level_var = level_var,
        serial = vars::serial.test(
            level_var, lags.pt = serial_lags, type = "PT.adjusted"
        ),
        arch = vars::arch.test(
            level_var, lags.multi = arch_lags, multivariate.only = TRUE
        ),
        roots = vars::roots(level_var, modulus = FALSE)
    )
}
