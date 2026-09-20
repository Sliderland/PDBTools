# Short diagnostic fit for the Ansley--Kohn VAR implementation.
#
# This script does not register anything in PosteriorDB and does not write
# reference-draw files. It uses the same small3 data construction as the
# interactive Heaps workflow, but runs only a short two-chain fit.

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
})

source(file.path("HeapsStanPrograms", "read.R"))

stan_file <- file.path("HeapsStanPrograms", "ansleykohn.stan")
output_dir <- file.path("ansleykohn_diagnostics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_prefix <- file.path(output_dir, "ansleykohn_small3_short")
fit_path <- paste0(output_prefix, "_fit.rds")
diagnostics_path <- paste0(output_prefix, "_diagnostics.rds")
report_path <- paste0(output_prefix, "_report.txt")
data_path <- paste0(output_prefix, "_data.rds")

# Match the established small3 preprocessing used by the batch workflow.
heaps_program_path <- normalizePath("HeapsStanPrograms", mustWork = TRUE)
y <- process_data(
  num = 3L,
  omit = c(1:2, 197:200),
  Nahead = 40L,
  data_dir_path = file.path(heaps_program_path, "data")
)$y

stan_data <- list(
  m = ncol(y),
  p = 4L,
  N = nrow(y),
  y = y
)

saveRDS(stan_data, data_path)

sampling_args <- list(
  chains = 2L,
  iter = 300L,
  warmup = 150L,
  refresh = 25L,
  seed = 20260920L,
  init = 0,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 12L
  ),
  cores = min(2L, parallel::detectCores())
)

warnings_seen <- character()
fit <- NULL
fit_error <- NULL

fit <- tryCatch(
  withCallingHandlers(
    rstan::sampling(
      object = rstan::stan_model(file = stan_file),
      data = stan_data,
      chains = sampling_args$chains,
      iter = sampling_args$iter,
      warmup = sampling_args$warmup,
      refresh = sampling_args$refresh,
      seed = sampling_args$seed,
      init = sampling_args$init,
      control = sampling_args$control,
      cores = sampling_args$cores,
      save_warmup = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) {
    fit_error <<- conditionMessage(e)
    NULL
  }
)

diagnostics <- list(
  timestamp = Sys.time(),
  stan_file = normalizePath(stan_file, mustWork = TRUE),
  data_dimensions = list(
    m = stan_data$m,
    p = stan_data$p,
    N = stan_data$N
  ),
  data_summary = list(
    finite = all(is.finite(stan_data$y)),
    range = range(stan_data$y),
    column_sd = apply(stan_data$y, 2L, sd)
  ),
  sampling_args = sampling_args,
  warnings = unique(warnings_seen),
  error = fit_error,
  fit_returned = !is.null(fit)
)

if (!is.null(fit)) {
  saveRDS(fit, fit_path)

  draws <- posterior::as_draws_array(fit)
  variable_names <- posterior::variables(draws)

  finite_by_variable <- vapply(seq_along(variable_names), function(i) {
    all(is.finite(as.numeric(draws[, , i])))
  }, logical(1))
  names(finite_by_variable) <- variable_names

  phi_names <- grep("^phi\\[", variable_names, value = TRUE)
  gamma_names <- grep("^Gamma\\[", variable_names, value = TRUE)

  sampler_params <- rstan::get_sampler_params(
    fit,
    inc_warmup = FALSE
  )

  diagnostics$draw_dimensions <- dim(draws)
  diagnostics$nonfinite_variables <- names(finite_by_variable)[
    !finite_by_variable
  ]
  diagnostics$phi_variables <- list(
    count = length(phi_names),
    all_finite = all(finite_by_variable[phi_names]),
    nonfinite = phi_names[!finite_by_variable[phi_names]]
  )
  diagnostics$Gamma_variables <- list(
    count = length(gamma_names),
    all_finite = all(finite_by_variable[gamma_names]),
    nonfinite = gamma_names[!finite_by_variable[gamma_names]]
  )
  diagnostics$stan_summary <- rstan::summary(fit)$summary
  diagnostics$sampler_summary <- lapply(sampler_params, function(x) {
    list(
      divergences = sum(x[, "divergent__"]),
      treedepth_max = max(x[, "treedepth__"]),
      treedepth_at_max = sum(x[, "treedepth__"] >= sampling_args$control$max_treedepth),
      energy_finite = all(is.finite(x[, "energy__"]))
    )
  })
}

saveRDS(diagnostics, diagnostics_path)

report_lines <- c(
  "Ansley--Kohn short diagnostic",
  paste("Timestamp:", diagnostics$timestamp),
  paste("Fit returned:", diagnostics$fit_returned),
  paste("Error:", if (is.null(diagnostics$error)) "<none>" else diagnostics$error),
  paste("Data finite:", diagnostics$data_summary$finite),
  paste("Data dimensions:", paste(unlist(diagnostics$data_dimensions), collapse = " x ")),
  paste("Warnings:", length(diagnostics$warnings))
)

if (!is.null(fit)) {
  report_lines <- c(
    report_lines,
    paste("Draw dimensions:", paste(diagnostics$draw_dimensions, collapse = " x ")),
    paste("Non-finite variables:", length(diagnostics$nonfinite_variables)),
    paste("All phi variables finite:", diagnostics$phi_variables$all_finite),
    paste("All Gamma variables finite:", diagnostics$Gamma_variables$all_finite),
    paste(
      "Total divergences:",
      sum(vapply(diagnostics$sampler_summary, `[[`, numeric(1), "divergences"))
    ),
    paste(
      "Maximum treedepth observed:",
      max(vapply(diagnostics$sampler_summary, `[[`, numeric(1), "treedepth_max"))
    )
  )
}

writeLines(report_lines, report_path)

message("Saved diagnostics to: ", diagnostics_path)
message("Saved report to: ", report_path)
if (!is.null(fit)) message("Saved short fit to: ", fit_path)
