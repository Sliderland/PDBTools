#!/usr/bin/env Rscript

# Launch the interactive 3D explorer for a saved rstan fit, including sampler
# diagnostics such as divergent transitions.
#
# Default fit:
#   Rscript visualize_failed_fit_with_divergences.R
#
# A different fit:
#   Rscript visualize_failed_fit_with_divergences.R \
#     --fit-path=/path/to/fit.rds

suppressPackageStartupMessages({
  library(posterior)
  library(rstan)
})

source("/Users/gerpr308/Documents/PDBTools/stan_3d_tools.R")

fit_path <- path.expand(
  paste0(
    "~/Documents/PDBTools_diagnostics/failed_reference_fits/",
    "heaps_med20-statprior_var_20260923_114608_pid95796.rds"
  )
)

args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (grepl("^--fit-path=", arg)) {
    fit_path <- path.expand(sub("^--fit-path=", "", arg))
  } else {
    stop("Unknown argument: ", arg, call. = FALSE)
  }
}

if (!file.exists(fit_path)) {
  stop("Fit file does not exist: ", fit_path, call. = FALSE)
}

fit <- readRDS(fit_path)
if (!inherits(fit, "stanfit")) {
  stop("The saved object is not an rstan stanfit object.", call. = FALSE)
}

diagnostics <- posterior::as_draws_df(
  rstan::get_sampler_params(fit, inc_warmup = FALSE)
)
divergent_count <- sum(diagnostics$divergent__ == 1, na.rm = TRUE)
transition_count <- sum(!is.na(diagnostics$divergent__))

cat("Fit:", fit_path, "\n")
cat("Post-warmup transitions:", transition_count, "\n")
cat("Divergent transitions:", divergent_count, "\n")

launch_stan_3d(fit)
