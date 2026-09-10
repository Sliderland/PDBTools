# Posterior prediction and residual diagnostics for Bayesian VAR(p) models.
#
# Array conventions:
# - phi_draws:  draws x lags x variables x variables
# - sigma_draws: draws x variables x variables
# - fitted/forecast draws: draws x time-or-horizon x variables
#
# The conditional mean matches the Heaps Stan programs:
#   mu + sum_l Phi_l %*% (y[t - l, ] - mu)
#
# Typical use with an rstan fit:
#   posterior <- extract_var_draws(fit)
#   fitted <- var_fitted_draws(posterior$phi_draws, y)
#   forecasts <- var_forecast_draws(
#       posterior$phi_draws, posterior$sigma_draws, y, horizon = 8L
#   )
#   checks <- posterior_residual_diagnostics(posterior$phi_draws, y)

validate_var_inputs <- function(phi_draws, y, sigma_draws = NULL, mu = NULL) {
    y <- as.matrix(y)
    phi_draws <- as.array(phi_draws)

    if (length(dim(phi_draws)) != 4L) {
        stop(
            "`phi_draws` must have dimensions draws x lags x variables x variables.",
            call. = FALSE
        )
    }
    if (!is.numeric(y) || !is.numeric(phi_draws)) {
        stop("`y` and `phi_draws` must be numeric.", call. = FALSE)
    }

    n_draws <- dim(phi_draws)[1L]
    p <- dim(phi_draws)[2L]
    k <- ncol(y)
    if (!identical(dim(phi_draws)[3:4], c(k, k))) {
        stop("The coefficient matrices do not match `ncol(y)`.", call. = FALSE)
    }
    if (nrow(y) <= p) {
        stop("`y` must contain more rows than the VAR lag order.", call. = FALSE)
    }
    if (anyNA(y) || anyNA(phi_draws) || any(!is.finite(y)) ||
        any(!is.finite(phi_draws))) {
        stop("`y` and `phi_draws` must contain only finite values.", call. = FALSE)
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
    }

    if (is.null(mu)) {
        mu <- numeric(k)
    }
    if (is.atomic(mu) && is.null(dim(mu)) && length(mu) == k) {
        mu <- matrix(rep(mu, times = n_draws), nrow = n_draws, byrow = TRUE)
    } else {
        mu <- as.matrix(mu)
    }
    if (!identical(dim(mu), c(n_draws, k))) {
        stop(
            "`mu` must be a length-k vector or a draws x k matrix.",
            call. = FALSE
        )
    }

    list(
        phi_draws = phi_draws,
        sigma_draws = sigma_draws,
        y = y,
        mu = mu,
        n_draws = n_draws,
        p = p,
        k = k
    )
}

select_posterior_draws <- function(n_draws, draws = NULL, max_draws = NULL) {
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

posterior_point_estimate <- function(draws, method = c("median", "mean")) {
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

extract_var_draws <- function(x, phi_name = "phi", sigma_name = "Sigma",
                              include_sigma = TRUE) {
    extracted <- if (inherits(x, "stanfit")) {
        if (!requireNamespace("rstan", quietly = TRUE)) {
            stop("The `rstan` package is required for a `stanfit` input.", call. = FALSE)
        }
        parameters <- c(phi_name, if (include_sigma) sigma_name)
        rstan::extract(x, pars = parameters, permuted = TRUE)
    } else if (is.list(x)) {
        x
    } else {
        stop("`x` must be a `stanfit` or an extracted-draw list.", call. = FALSE)
    }

    if (is.null(extracted[[phi_name]])) {
        stop("No coefficient draws named `", phi_name, "` were found.", call. = FALSE)
    }
    result <- list(phi_draws = extracted[[phi_name]])
    if (include_sigma) {
        if (is.null(extracted[[sigma_name]])) {
            stop("No covariance draws named `", sigma_name, "` were found.", call. = FALSE)
        }
        result$sigma_draws <- extracted[[sigma_name]]
    }
    result
}

var_fitted_draws <- function(phi_draws, y, mu = NULL, draws = NULL,
                             max_draws = NULL) {
    inputs <- validate_var_inputs(phi_draws, y, mu = mu)
    selected <- select_posterior_draws(
        inputs$n_draws,
        draws = draws,
        max_draws = max_draws
    )

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
        for (t in seq.int(inputs$p + 1L, nrow(inputs$y))) {
            conditional_mean <- inputs$mu[d, ]
            for (lag in seq_len(inputs$p)) {
                conditional_mean <- conditional_mean +
                    inputs$phi_draws[d, lag, , ] %*%
                        (inputs$y[t - lag, ] - inputs$mu[d, ])
            }
            fitted[d_out, t - inputs$p, ] <- conditional_mean
        }
    }
    fitted
}

var_fitted_point <- function(phi_draws, y, mu = NULL,
                             method = c("median", "mean")) {
    method <- match.arg(method)
    phi <- posterior_point_estimate(phi_draws, method = method)
    phi <- array(phi, dim = c(1L, dim(phi)))
    point_mu <- if (is.null(mu) || is.atomic(mu) && is.null(dim(mu))) {
        mu
    } else {
        posterior_point_estimate(mu, method = method)
    }
    drop(var_fitted_draws(phi, y, mu = point_mu)[1L, , , drop = FALSE])
}

var_residual_draws <- function(phi_draws, y, mu = NULL, draws = NULL,
                               max_draws = NULL) {
    fitted <- var_fitted_draws(
        phi_draws,
        y,
        mu = mu,
        draws = draws,
        max_draws = max_draws
    )
    p <- dim(as.array(phi_draws))[2L]
    sweep(fitted, c(2L, 3L), as.matrix(y)[-(seq_len(p)), , drop = FALSE], "-") * -1
}

rmvn_chol <- function(sigma) {
    drop(t(chol(sigma)) %*% stats::rnorm(nrow(sigma)))
}

var_forecast_draws <- function(phi_draws, sigma_draws, y, horizon,
                               mu = NULL, draws = NULL, max_draws = NULL,
                               seed = NULL) {
    inputs <- validate_var_inputs(phi_draws, y, sigma_draws, mu)
    if (length(horizon) != 1L || !is.finite(horizon) || horizon < 1L ||
        horizon %% 1 != 0) {
        stop("`horizon` must be a positive integer.", call. = FALSE)
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    selected <- select_posterior_draws(
        inputs$n_draws,
        draws = draws,
        max_draws = max_draws
    )

    forecasts <- array(
        NA_real_,
        dim = c(length(selected), as.integer(horizon), inputs$k),
        dimnames = list(
            draw = selected,
            horizon = seq_len(horizon),
            variable = colnames(inputs$y)
        )
    )

    for (d_out in seq_along(selected)) {
        d <- selected[[d_out]]
        history <- inputs$y
        for (step in seq_len(horizon)) {
            conditional_mean <- inputs$mu[d, ]
            for (lag in seq_len(inputs$p)) {
                conditional_mean <- conditional_mean +
                    inputs$phi_draws[d, lag, , ] %*%
                        (history[nrow(history) - lag + 1L, ] - inputs$mu[d, ])
            }
            next_value <- drop(conditional_mean) +
                rmvn_chol(inputs$sigma_draws[d, , ])
            forecasts[d_out, step, ] <- next_value
            history <- rbind(history, next_value)
        }
    }
    forecasts
}

summarize_posterior_array <- function(x, probs = c(0.05, 0.5, 0.95)) {
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

residual_acf <- function(residuals, lag_max = 20L,
                         correlation = TRUE) {
    residuals <- as.matrix(residuals)
    if (nrow(residuals) <= lag_max) {
        stop("`lag_max` must be smaller than the number of residuals.", call. = FALSE)
    }

    marginal <- lapply(seq_len(ncol(residuals)), function(j) {
        drop(stats::acf(
            residuals[, j],
            lag.max = lag_max,
            plot = FALSE,
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

multivariate_portmanteau <- function(residuals, lags = 16L,
                                     fitted_lags = 0L,
                                     adjusted = TRUE) {
    input_name <- deparse(substitute(residuals))
    residuals <- as.matrix(residuals)
    n <- nrow(residuals)
    k <- ncol(residuals)
    if (lags <= fitted_lags || lags >= n) {
        stop("Require `fitted_lags < lags < nrow(residuals)`.", call. = FALSE)
    }

    c0 <- crossprod(residuals) / n
    c0_inverse <- tryCatch(
        solve(c0),
        error = function(error) {
            stop("Residual covariance is singular.", call. = FALSE)
        }
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
    degrees_freedom <- k^2 * (lags - fitted_lags)
    structure(
        list(
            statistic = c("Chi-squared" = statistic),
            parameter = c(df = degrees_freedom),
            p.value = stats::pchisq(
                statistic,
                df = degrees_freedom,
                lower.tail = FALSE
            ),
            method = paste(
                if (adjusted) "Adjusted" else "Asymptotic",
                "multivariate Portmanteau test"
            ),
            data.name = input_name
        ),
        class = "htest"
    )
}

vech_rows <- function(residuals) {
    residuals <- as.matrix(residuals)
    keep <- lower.tri(matrix(0, ncol(residuals), ncol(residuals)), diag = TRUE)
    t(vapply(seq_len(nrow(residuals)), function(i) {
        tcrossprod(residuals[i, ])[keep]
    }, numeric(sum(keep))))
}

multivariate_arch_lm <- function(residuals, lags = 1L) {
    input_name <- deparse(substitute(residuals))
    residuals <- scale(as.matrix(residuals))
    n <- nrow(residuals)
    k <- ncol(residuals)
    if (lags < 1L || lags >= n) {
        stop("Require `1 <= lags < nrow(residuals)`.", call. = FALSE)
    }

    squared_products <- vech_rows(residuals)
    embedded <- embed(squared_products, lags + 1L)
    n_products <- ncol(squared_products)
    response <- embedded[, seq_len(n_products), drop = FALSE]
    predictors <- embedded[, -seq_len(n_products), drop = FALSE]
    if (nrow(response) <= ncol(predictors) + 1L) {
        stop(
            paste0(
                "The multivariate ARCH regression is underdetermined: ",
                nrow(response), " rows for ", ncol(predictors) + 1L,
                " coefficients per response. Reduce `lags` or use more data."
            ),
            call. = FALSE
        )
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
            p.value = stats::pchisq(
                statistic,
                df = degrees_freedom,
                lower.tail = FALSE
            ),
            r.squared = r_squared,
            method = "Multivariate ARCH-LM test",
            data.name = input_name
        ),
        class = "htest"
    )
}

posterior_residual_diagnostics <- function(
    phi_draws,
    y,
    mu = NULL,
    draws = NULL,
    max_draws = 200L,
    serial_lags = 16L,
    arch_lags = 1L,
    alpha = 0.05
) {
    residual_draws <- var_residual_draws(
        phi_draws,
        y,
        mu = mu,
        draws = draws,
        max_draws = max_draws
    )
    fitted_lags <- dim(as.array(phi_draws))[2L]

    diagnostics <- lapply(seq_len(dim(residual_draws)[1L]), function(d) {
        serial <- multivariate_portmanteau(
            residual_draws[d, , ],
            lags = serial_lags,
            fitted_lags = fitted_lags,
            adjusted = TRUE
        )
        arch <- multivariate_arch_lm(
            residual_draws[d, , ],
            lags = arch_lags
        )
        list(serial = serial, arch = arch)
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
            proportion_rejecting = c(mean(serial_p < alpha), mean(arch_p < alpha)),
            alpha = alpha,
            row.names = NULL
        )
    )
}

classical_var_diagnostics <- function(y, p, type = "none",
                                      serial_lags = 16L,
                                      arch_lags = 5L) {
    if (!requireNamespace("vars", quietly = TRUE)) {
        stop(
            "Install the `vars` package to use `classical_var_diagnostics()`.",
            call. = FALSE
        )
    }
    fit <- vars::VAR(as.matrix(y), p = p, type = type)
    list(
        fit = fit,
        serial = vars::serial.test(
            fit,
            lags.pt = serial_lags,
            type = "PT.adjusted"
        ),
        arch = vars::arch.test(
            fit,
            lags.multi = arch_lags,
            multivariate.only = TRUE
        )
    )
}
