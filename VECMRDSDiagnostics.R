# Batch diagnostics for the VECM_info.RDS archive produced by
# unconstrainedVECM.R.
#
# This file orchestrates VECMPosteriorDiagnostics.R. It reconstructs the known
# urca::denmark input for vecm_urca_hmc and accepts additional Stan data lists
# for the simulated-data fits when they become available.
#
# For the recovered archive in this project:
#   source("VECMRDSDiagnostics.R")
#   recovered <- readRDS("recovered_vecm_data.rds")
#   checks <- diagnose_vecm_rds(
#       stan_data_by_model = recovered$stan_data,
#       max_draws = 200L
#   )
#   summary(checks)
#   checks$vecm_urca_hmc$residual_summary
#   checks$vecm_urca_hmc$innovation_ppc
#
# The synthetic data in that file omit the unrecovered final observation.
# Predictive calibration below is in-sample and conditional on observed lags;
# it is not a held-out or leave-one-out forecast score.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Source relative to this file when possible, otherwise relative to the current
# working directory.
local({
    caller_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    source_path <- if (is.null(caller_file)) {
        "VECMPosteriorDiagnostics.R"
    } else {
        file.path(dirname(normalizePath(caller_file)), "VECMPosteriorDiagnostics.R")
    }
    if (!exists("extract_saved_vecm_draws", mode = "function")) {
        source(source_path)
    }
})

default_vecm_archive_path <- function() {
    file.path(
        "/Users/gerpr308/Library/CloudStorage",
        "OneDrive-Uppsalauniversitet",
        "Bayesian Time Series",
        "Bayesian Time Series Code",
        "VECM_info.RDS"
    )
}

reconstruct_urca_denmark_stan_data <- function(p = 4L, rank = 2L,
                                               seasonality = 4L) {
    if (!requireNamespace("urca", quietly = TRUE)) {
        stop("Install the `urca` package to reconstruct the Denmark data.", call. = FALSE)
    }
    environment <- new.env(parent = emptyenv())
    utils::data("denmark", package = "urca", envir = environment)
    if (!exists("denmark", envir = environment, inherits = FALSE)) {
        stop("Could not load `urca::denmark`.", call. = FALSE)
    }
    denmark <- get("denmark", envir = environment)
    y <- as.matrix(denmark[, -1L, drop = FALSE])
    seasons <- factor(rep(seq_len(seasonality), length.out = nrow(y)))
    raw_dummies <- stats::model.matrix(~ seasons - 1L)
    centered <- scale(raw_dummies, center = TRUE, scale = FALSE)
    seasonal_dummies <- centered[, seq_len(seasonality - 1L), drop = FALSE]
    storage.mode(seasonal_dummies) <- "double"

    list(
        t = nrow(y),
        N = ncol(y),
        p = as.integer(p),
        h = as.integer(rank),
        y = y,
        num_exo = ncol(seasonal_dummies),
        D = seasonal_dummies,
        error_correction_lag = as.integer(p),
        specification = "long_run_with_exogenous"
    )
}

default_vecm_stan_data <- function() {
    list(vecm_urca_hmc = reconstruct_urca_denmark_stan_data())
}

validate_vecm_stan_data <- function(stan_data, model_name) {
    if (!is.list(stan_data) || is.null(stan_data$y) || is.null(stan_data$p)) {
        stop(
            "Stan data for `", model_name, "` must contain at least `y` and `p`.",
            call. = FALSE
        )
    }
    stan_data$y <- as.matrix(stan_data$y)
    stan_data$p <- as.integer(stan_data$p)
    if (nrow(stan_data$y) <= stan_data$p) {
        stop("The supplied data contain too few observations.", call. = FALSE)
    }
    if (is.null(stan_data$error_correction_lag)) {
        stan_data$error_correction_lag <- if (
            model_name %in% c("vecm_long_run", "vecm_urca_hmc")
        ) stan_data$p else 1L
    }
    stan_data$error_correction_lag <- as.integer(stan_data$error_correction_lag)
    if (stan_data$error_correction_lag < 1L ||
        stan_data$error_correction_lag > stan_data$p) {
        stop("`error_correction_lag` must lie between 1 and `p`.", call. = FALSE)
    }
    if (!is.null(stan_data$D)) {
        stan_data$D <- as.matrix(stan_data$D)
        if (nrow(stan_data$D) != nrow(stan_data$y)) {
            stop("`D` and `y` must have the same number of rows.", call. = FALSE)
        }
    }
    stan_data
}

saved_vecm_fitted_draws <- function(parameters, stan_data, draws = NULL,
                                    max_draws = 200L) {
    stan_data <- validate_vecm_stan_data(stan_data, parameters$model_name)
    inputs <- validate_vecm_inputs(
        parameters$alpha_draws,
        parameters$beta_draws,
        parameters$gamma_draws,
        stan_data$y,
        intercept = parameters$intercept_draws
    )
    if (inputs$p != stan_data$p) {
        stop("The draw lag order does not match `stan_data$p`.", call. = FALSE)
    }
    selected <- select_vecm_draws(inputs$n_draws, draws, max_draws)
    exogenous <- parameters$exogenous_draws
    if (!is.null(stan_data$D) && is.null(exogenous)) {
        stop("The data include `D`, but no saved `phi` draws were found.", call. = FALSE)
    }
    if (!is.null(exogenous)) {
        exogenous <- as.array(exogenous)
        if (is.null(stan_data$D)) {
            stop("Saved `phi` draws require an exogenous matrix `D`.", call. = FALSE)
        }
        if (!identical(
            dim(exogenous),
            c(inputs$n_draws, inputs$k, ncol(stan_data$D))
        )) {
            stop("The saved `phi` draws do not match `D`.", call. = FALSE)
        }
    }

    y <- stan_data$y
    dy <- diff(y)
    fitted <- array(
        NA_real_,
        dim = c(length(selected), nrow(y) - inputs$p, inputs$k),
        dimnames = list(
            draw = selected,
            time = seq.int(inputs$p + 1L, nrow(y)),
            variable = colnames(y)
        )
    )
    for (d_out in seq_along(selected)) {
        d <- selected[[d_out]]
        alpha <- matrix(inputs$alpha_draws[d, , ], nrow = inputs$k)
        beta <- matrix(inputs$beta_draws[d, , ], nrow = inputs$k)
        pi <- alpha %*% t(beta)
        for (t in seq.int(inputs$p + 1L, nrow(y))) {
            mean_change <- inputs$intercept[d, ] +
                pi %*% y[t - stan_data$error_correction_lag, ]
            for (lag in seq_len(inputs$q)) {
                mean_change <- mean_change +
                    inputs$gamma_draws[d, lag, , ] %*% dy[t - lag - 1L, ]
            }
            if (!is.null(exogenous)) {
                mean_change <- mean_change +
                    exogenous[d, , ] %*% stan_data$D[t, ]
            }
            fitted[d_out, t - inputs$p, ] <- mean_change
        }
    }
    fitted
}

saved_vecm_residual_draws <- function(parameters, stan_data, draws = NULL,
                                      max_draws = 200L) {
    fitted <- saved_vecm_fitted_draws(parameters, stan_data, draws, max_draws)
    p <- as.integer(stan_data$p)
    observed <- diff(as.matrix(stan_data$y))[-seq_len(p - 1L), , drop = FALSE]
    sweep(fitted, c(2L, 3L), observed, "-") * -1
}

summarize_saved_vecm_sampler <- function(archive, model_name) {
    diagnostics <- archive$diag[[model_name]]
    draws <- archive$samps[[model_name]]
    variables <- posterior::variables(draws)
    primary <- c("alpha", "beta", "mu", "xi", "gamma", "phi",
                 "L", "L_Omega", "L_sigma")
    selected <- variables[vapply(variables, function(variable) {
        any(startsWith(variable, paste0(primary, "[")))
    }, logical(1))]
    if (!length(selected)) {
        stop("No primary VECM parameters were found in the saved draws.",
             call. = FALSE)
    }
    summary <- posterior::summarise_draws(
        posterior::subset_draws(draws, variable = selected),
        rhat = posterior::rhat,
        ess_bulk = posterior::ess_bulk,
        ess_tail = posterior::ess_tail
    )
    list(
        divergences = sum(diagnostics$num_divergent),
        max_treedepth_hits = sum(diagnostics$num_max_treedepth),
        minimum_ebfmi = min(diagnostics$ebfmi),
        maximum_rhat = max(summary$rhat, na.rm = TRUE),
        minimum_bulk_ess = min(summary$ess_bulk, na.rm = TRUE),
        minimum_tail_ess = min(summary$ess_tail, na.rm = TRUE)
    )
}

# Conditional, in-sample posterior predictive check. The observed innovation
# statistic is compared with a fresh Gaussian innovation series for each draw.
# This is not a held-out or leave-one-out predictive assessment.
saved_vecm_innovation_ppc <- function(residual_draws, sigma_draws,
                                      lag_max = 16L, seed = 123L) {
    n_draws <- dim(residual_draws)[1L]
    n <- dim(residual_draws)[2L]
    k <- dim(residual_draws)[3L]
    selected <- as.integer(dimnames(residual_draws)[[1L]])
    if (length(selected) != n_draws || anyNA(selected)) {
        stop("Residual draws must retain their original draw indices.",
             call. = FALSE)
    }
    if (lag_max < 1L || lag_max >= n) {
        stop("`lag_max` must lie between 1 and nrow(residuals) - 1.",
             call. = FALSE)
    }
    if (!is.null(seed)) {
        had_seed <- exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)
        if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit({
            if (had_seed) {
                assign(".Random.seed", previous_seed, envir = .GlobalEnv)
            } else if (exists(".Random.seed", envir = .GlobalEnv,
                              inherits = FALSE)) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        })
        set.seed(seed)
    }
    max_acf <- function(x) {
        max(vapply(seq_len(ncol(x)), function(j) {
            max(abs(stats::acf(x[, j], lag.max = lag_max,
                               plot = FALSE)$acf[-1L]))
        }, numeric(1)))
    }
    max_squared_lag1_acf <- function(x) {
        max(vapply(seq_len(ncol(x)), function(j) {
            stats::acf(x[, j]^2, lag.max = 1L, plot = FALSE)$acf[2L]
        }, numeric(1)))
    }
    observed <- replicate <- matrix(NA_real_, nrow = n_draws, ncol = 2L)
    for (d in seq_len(n_draws)) {
        sigma <- sigma_draws[selected[d], , ]
        simulated <- matrix(stats::rnorm(n * k), nrow = n, ncol = k) %*%
            chol(sigma)
        residuals <- matrix(residual_draws[d, , ], nrow = n, ncol = k)
        observed[d, ] <- c(max_acf(residuals),
                           max_squared_lag1_acf(residuals))
        replicate[d, ] <- c(max_acf(simulated),
                            max_squared_lag1_acf(simulated))
    }
    data.frame(
        check = c("maximum_absolute_residual_acf",
                  "maximum_squared_residual_lag1_acf"),
        observed_median = apply(observed, 2L, stats::median),
        replicated_median = apply(replicate, 2L, stats::median),
        posterior_predictive_tail_probability = colMeans(replicate >= observed),
        row.names = NULL
    )
}

saved_vecm_predictive_calibration <- function(residual_draws, sigma_draws) {
    selected <- as.integer(dimnames(residual_draws)[[1L]])
    k <- dim(residual_draws)[3L]
    names <- dimnames(residual_draws)[[3L]]
    if (is.null(names)) names <- paste0("variable_", seq_len(k))
    data.frame(
        variable = names,
        central_90_coverage = vapply(seq_len(k), function(j) {
            standardized <- sweep(
                matrix(residual_draws[, , j], nrow = length(selected)), 1L,
                sqrt(sigma_draws[selected, j, j]), "/"
            )
            pit <- colMeans(stats::pnorm(standardized))
            mean(pit > 0.05 & pit < 0.95)
        }, numeric(1)),
        row.names = NULL
    )
}

diagnose_saved_vecm_model <- function(archive, model_name, stan_data = NULL,
                                      max_draws = 200L, serial_lags = 16L,
                                      arch_lags = 1L,
                                      rank_tolerance = 1e-8,
                                      unit_root_tolerance = 1e-6) {
    parameters <- suppressWarnings(extract_saved_vecm_draws(archive, model_name))
    pi_draws <- vecm_pi_draws(parameters$alpha_draws, parameters$beta_draws)
    pi_ranks <- vapply(seq_len(dim(pi_draws)[1L]), function(d) {
        values <- svd(pi_draws[d, , ], nu = 0L, nv = 0L)$d
        sum(values > rank_tolerance * max(1, values[[1L]]))
    }, integer(1))
    result <- list(
        model_name = model_name,
        sampler_diagnostics = archive$diag[[model_name]],
        sampler_summary = summarize_saved_vecm_sampler(archive, model_name),
        nuts_parameters = archive$np[[model_name]],
        parameters = parameters,
        parameter_checks = list(
            expected_rank = dim(parameters$alpha_draws)[3L],
            numerical_rank = pi_ranks,
            proportion_rank_ok = mean(pi_ranks == dim(parameters$alpha_draws)[3L]),
            median_pi = apply(pi_draws, c(2L, 3L), stats::median),
            median_cointegration_projection = apply(
                vecm_cointegration_projection_draws(parameters$beta_draws),
                c(2L, 3L), stats::median
            )
        ),
        status = "parameter diagnostics completed"
    )

    is_transitory <- !model_name %in% c("vecm_long_run", "vecm_urca_hmc")
    if (is_transitory) {
        result$rank_root_checks <- vecm_rank_root_diagnostics(
            parameters$alpha_draws,
            parameters$beta_draws,
            parameters$gamma_draws,
            rank_tolerance = rank_tolerance,
            unit_root_tolerance = unit_root_tolerance
        )
    }
    if (is.null(stan_data)) {
        result$status <- paste(
            result$status,
            "only; exact Stan data were not supplied"
        )
        return(result)
    }

    stan_data <- validate_vecm_stan_data(stan_data, model_name)
    residuals <- saved_vecm_residual_draws(
        parameters, stan_data, max_draws = max_draws
    )
    q <- dim(parameters$gamma_draws)[2L]
    per_draw <- lapply(seq_len(dim(residuals)[1L]), function(d) {
        list(
            serial = vecm_multivariate_portmanteau(
                residuals[d, , ], lags = serial_lags,
                fitted_difference_lags = q
            ),
            arch = vecm_multivariate_arch_lm(residuals[d, , ], lags = arch_lags)
        )
    })
    serial_p <- vapply(per_draw, function(x) x$serial$p.value, numeric(1))
    arch_p <- vapply(per_draw, function(x) x$arch$p.value, numeric(1))
    result$residual_draws <- residuals
    result$residual_diagnostics <- per_draw
    result$residual_summary <- data.frame(
        test = c("serial_correlation", "multivariate_arch"),
        median_p_value = c(stats::median(serial_p), stats::median(arch_p)),
        proportion_rejecting_0.05 = c(mean(serial_p < 0.05), mean(arch_p < 0.05)),
        row.names = NULL
    )
    result$innovation_ppc <- saved_vecm_innovation_ppc(
        residuals, parameters$sigma_draws,
        lag_max = min(serial_lags, dim(residuals)[2L] - 1L)
    )
    result$predictive_calibration <- saved_vecm_predictive_calibration(
        residuals, parameters$sigma_draws
    )
    result$status <- "parameter and residual diagnostics completed"
    result
}

diagnose_vecm_rds <- function(
    archive_path = default_vecm_archive_path(),
    models = NULL,
    stan_data_by_model = default_vecm_stan_data(),
    max_draws = 200L,
    serial_lags = 16L,
    arch_lags = 1L,
    keep_parameters = FALSE
) {
    archive <- if (is.character(archive_path)) readRDS(archive_path) else archive_path
    available <- intersect(names(archive$models), names(archive$samps))
    if (is.null(models)) {
        models <- available
    }
    unknown <- setdiff(models, available)
    if (length(unknown)) {
        stop("Unknown models: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }

    results <- setNames(vector("list", length(models)), models)
    for (model_name in models) {
        message("Diagnosing ", model_name, "...")
        results[[model_name]] <- tryCatch(
            diagnose_saved_vecm_model(
                archive = archive,
                model_name = model_name,
                stan_data = stan_data_by_model[[model_name]],
                max_draws = max_draws,
                serial_lags = serial_lags,
                arch_lags = arch_lags
            ),
            error = function(error) list(
                model_name = model_name,
                status = "failed",
                error = conditionMessage(error)
            )
        )
        if (!keep_parameters) {
            results[[model_name]]$parameters <- NULL
            results[[model_name]]$residual_draws <- NULL
            results[[model_name]]$nuts_parameters <- NULL
        }
    }
    class(results) <- c("vecm_rds_diagnostics", "list")
    results
}

summary.vecm_rds_diagnostics <- function(object, ...) {
    data.frame(
        model = names(object),
        status = vapply(object, function(x) x$status %||% "unknown", character(1)),
        error = vapply(object, function(x) x$error %||% "", character(1)),
        rank_ok = vapply(object, function(x) {
            x$parameter_checks$proportion_rank_ok %||% NA_real_
        }, numeric(1)),
        divergences = vapply(object, function(x) {
            x$sampler_summary$divergences %||% NA_real_
        }, numeric(1)),
        treedepth_hits = vapply(object, function(x) {
            x$sampler_summary$max_treedepth_hits %||% NA_real_
        }, numeric(1)),
        max_rhat = vapply(object, function(x) {
            x$sampler_summary$maximum_rhat %||% NA_real_
        }, numeric(1)),
        min_ebfmi = vapply(object, function(x) {
            x$sampler_summary$minimum_ebfmi %||% NA_real_
        }, numeric(1)),
        min_bulk_ess = vapply(object, function(x) {
            x$sampler_summary$minimum_bulk_ess %||% NA_real_
        }, numeric(1)),
        residual_checks = vapply(object, function(x) {
            !is.null(x$residual_summary)
        }, logical(1)),
        row.names = NULL
    )
}
