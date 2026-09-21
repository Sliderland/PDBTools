library(posterior)

archive_path <- paste0(
  "/Users/gerpr308/Library/CloudStorage/",
  "OneDrive-Uppsalauniversitet/Bayesian Time Series/",
  "Bayesian Time Series Code/VECM_info.RDS"
)
archive <- readRDS(archive_path) # Approximately 1.7 GB

# Use separated draws from one chain. Increase this if desired.
iterations <- unique(as.integer(round(seq(1, 1000, length.out = 80))))

get_arrays <- function(model, families) {
  draws <- archive$samps[[model]]
  variable_names <- dimnames(draws)[[3]]

  columns <- which(vapply(
    variable_names,
    function(v) any(startsWith(v, paste0(families, "["))),
    logical(1)
  ))

  selected <- posterior::as_draws_array(
    draws[iterations, 1L, columns, drop = FALSE]
  )
  lapply(posterior::as_draws_rvars(selected), posterior::draws_of)
}

recover_ect_levels <- function(z) {
  n_draws <- dim(z$pi)[1]
  k <- dim(z$pi)[2]

  # For every draw: t(pi) %*% y_t = ect_t
  A <- do.call(
    rbind,
    lapply(
      seq_len(n_draws),
      function(d) t(z$pi[d, , ])
    )
  )
  B <- do.call(
    rbind,
    lapply(
      seq_len(n_draws),
      function(d) t(z$ect[d, , ])
    )
  )

  if (qr(A)$rank != k) {
    stop("The stacked ECT system lacks full rank.")
  }

  y <- t(qr.solve(A, B))
  list(
    y = y,
    rank = qr(A)$rank,
    max_equation_error = max(abs(A %*% t(y) - B))
  )
}

recover_lagged_differences <- function(z) {
  n_draws <- dim(z$topblock)[1]
  n_unknowns <- dim(z$topblock)[3]

  # For every draw: topblock %*% dX_t = short_run_t
  A <- do.call(
    rbind,
    lapply(
      seq_len(n_draws),
      function(d) z$topblock[d, , ]
    )
  )
  B <- do.call(
    rbind,
    lapply(
      seq_len(n_draws),
      function(d) t(z$short_run[d, , ])
    )
  )

  if (qr(A)$rank != n_unknowns) {
    stop("The stacked short-run system lacks full rank.")
  }

  dX <- t(qr.solve(A, B))
  list(
    dX = dX,
    rank = qr(A)$rank,
    max_equation_error = max(abs(A %*% t(dX) - B))
  )
}

baseline_arrays <- get_arrays(
  "vecm_baseline",
  c("pi", "ect", "topblock", "short_run")
)
priors_arrays <- get_arrays("vecm_priors", c("pi", "ect"))
long_run_arrays <- get_arrays("vecm_long_run", c("pi", "ect"))

baseline <- recover_ect_levels(baseline_arrays) # y[4:999, ]
priors <- recover_ect_levels(priors_arrays) # y[4:999, ]
long_run <- recover_ect_levels(long_run_arrays) # y[1:996, ]
short_run <- recover_lagged_differences(baseline_arrays)

synthetic_y <- matrix(NA_real_, nrow = 1000, ncol = 5)
synthetic_y[1:996, ] <- long_run$y
synthetic_y[4:999, ] <- baseline$y
colnames(synthetic_y) <- c("LRM", "LRY", "LPY", "IBO", "IDE")

# Check the recovered levels against independently recovered dX.
implied_dX <- cbind(
  synthetic_y[4:999, ] - synthetic_y[3:998, ],
  synthetic_y[3:998, ] - synthetic_y[2:997, ],
  synthetic_y[2:997, ] - synthetic_y[1:996, ]
)

recovery_checks <- c(
  baseline_rank = baseline$rank,
  priors_rank = priors$rank,
  long_run_rank = long_run$rank,
  dX_rank = short_run$rank,
  baseline_equation_error = baseline$max_equation_error,
  priors_equation_error = priors$max_equation_error,
  long_run_equation_error = long_run$max_equation_error,
  dX_equation_error = short_run$max_equation_error,
  baseline_vs_priors = max(abs(baseline$y - priors$y)),
  baseline_vs_long_run = max(
    abs(baseline$y[1:993, ] - long_run$y[4:996, ])
  ),
  difference_consistency = max(abs(implied_dX - short_run$dX))
)
print(recovery_checks)

# Recreate the original observed-data input for vecm_urca_hmc.
denmark_env <- new.env(parent = emptyenv())
utils::data("denmark", package = "urca", envir = denmark_env)
denmark_y <- as.matrix(denmark_env$denmark[, -1L, drop = FALSE])

seasons <- factor(rep(1:4, length.out = nrow(denmark_y)))
seasonal_dummies <- scale(
  stats::model.matrix(~ seasons - 1),
  center = TRUE,
  scale = FALSE
)[, 1:3, drop = FALSE]

# All four reconstructed input datasets. The final synthetic row remains NA.
recovered_data <- list(
  vecm_baseline = synthetic_y,
  vecm_priors = synthetic_y,
  vecm_long_run = synthetic_y,
  vecm_urca_hmc = denmark_y
)

# For residual checks, use only the recovered 999 synthetic rows.
stan_data_by_model <- list(
  vecm_baseline = list(
    y = synthetic_y[1:999, , drop = FALSE],
    p = 4L,
    h = 2L,
    error_correction_lag = 1L
  ),
  vecm_priors = list(
    y = synthetic_y[1:999, , drop = FALSE],
    p = 4L,
    h = 2L,
    error_correction_lag = 1L
  ),
  vecm_long_run = list(
    y = synthetic_y[1:999, , drop = FALSE],
    p = 4L,
    h = 2L,
    error_correction_lag = 4L
  ),
  vecm_urca_hmc = list(
    y = denmark_y,
    p = 4L,
    h = 2L,
    D = seasonal_dummies,
    num_exo = 3L,
    error_correction_lag = 4L
  )
)

# Optional: preserve the result without changing the original archive.
saveRDS(
  list(
    data = recovered_data,
    stan_data = stan_data_by_model,
    checks = recovery_checks
  ),
  "/Users/gerpr308/Documents/PDBTools/recovered_vecm_data.rds"
)


source("/Users/gerpr308/Documents/PDBTools/VECMRDSDiagnostics.R")

recovered <- readRDS(
  "/Users/gerpr308/Documents/PDBTools/recovered_vecm_data.rds"
)

results <- diagnose_vecm_rds(
  stan_data_by_model = recovered$stan_data,
  max_draws = 200L,
  serial_lags = 16L,
  arch_lags = 1L
)

summary(results)

lapply(results, function(x) x$residual_summary)

lapply(results, function(x) {
  if (is.null(x$rank_root_checks)) {
    return(NULL)
  }
  colMeans(x$rank_root_checks[
    c("rank_ok", "unit_roots_ok", "no_explosive_roots")
  ])
})
source("/Users/gerpr308/Documents/PDBTools/VECMRDSDiagnostics.R")

recovered <- readRDS(
  "/Users/gerpr308/Documents/PDBTools/recovered_vecm_data.rds"
)

checks <- diagnose_vecm_rds(
  stan_data_by_model = recovered$stan_data,
  max_draws = 200L,
  serial_lags = 16L,
  arch_lags = 1L
)

summary(checks)

lapply(checks, `[[`, "sampler_summary")
lapply(checks, `[[`, "residual_summary")
lapply(checks, `[[`, "innovation_ppc")
lapply(checks, `[[`, "predictive_calibration")
source("/Users/gerpr308/Documents/PDBTools/VECMRDSDiagnostics.R")

recovered <- readRDS(
  "/Users/gerpr308/Documents/PDBTools/recovered_vecm_data.rds"
)

checks <- diagnose_vecm_rds(
  stan_data_by_model = recovered$stan_data,
  max_draws = 200L
)

summary(checks)
checks$vecm_urca_hmc$residual_summary
checks$vecm_urca_hmc$innovation_ppc
checks$vecm_urca_hmc$predictive_calibration
