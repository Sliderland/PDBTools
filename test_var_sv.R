# Sampling and diagnostics workflow for VAR_SV.stan.
# Run with: Rscript test_var_sv.R

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

stan_file <- "/Users/gerpr308/Library/CloudStorage/OneDrive-Uppsalauniversitet/Bayesian Time Series/Bayesian Time Series Code/VAR_SV.stan"
data_file <- "/Users/gerpr308/Library/CloudStorage/OneDrive-Uppsalauniversitet/Bayesian Time Series/Bayesian Time Series Code/sw2001.txt"
# Keep this retry separate from the previous fit.
output_dir <- file.path(getwd(), "var_sv_retry_095")
p <- 4L
chains <- 4L
parallel_chains <- 4L
iter_warmup <- 1000L
iter_sampling <- 20000L
thin <- 2L
adapt_delta <- 0.95
max_treedepth <- 12L
seed <- 20260917L

if (!file.exists(stan_file)) stop("Missing Stan file: ", stan_file)
if (!file.exists(data_file)) stop("Missing data file: ", data_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.delim(data_file, check.names = FALSE)
if (ncol(raw) < 2L) stop("Expected a date column and numeric series columns")
y <- as.matrix(raw[, -1L, drop = FALSE])
storage.mode(y) <- "double"
if (any(!is.finite(y))) stop("Data contain missing or non-finite values")
if (nrow(y) <= p) stop("The number of observations must exceed p")

stan_data <- list(
  t = nrow(y),
  N = ncol(y),
  p = p,
  y = y
)

message("Compiling VAR_SV.stan...")
model <- cmdstan_model(stan_file)

message("Sampling VAR-SV model...")
fit <- model$sample(
  data = stan_data,
  seed = seed,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  thin = thin,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  refresh = 100,
  output_dir = output_dir,
  show_messages = TRUE
)

sampler_draws <- fit$sampler_diagnostics(format = "draws_array")
sampler <- as_draws_df(sampler_draws)
divergences <- sum(sampler$divergent__ > 0)
treedepth_hits <- sum(sampler$treedepth__ >= max_treedepth)

ebfmi <- vapply(split(sampler$energy__, sampler$.chain), function(x) {
  if (length(x) < 2L || var(x) == 0) return(NA_real_)
  mean(diff(x)^2) / var(x)
}, numeric(1))
names(ebfmi) <- paste0("chain_", names(ebfmi))

draws <- fit$draws(inc_warmup = FALSE)
parameter_summary <- summarise_draws(
  draws,
  mean,
  sd,
  rhat = posterior::rhat,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail
)
bad_rhat <- parameter_summary[is.finite(parameter_summary$rhat) & parameter_summary$rhat > 1.01, , drop = FALSE]
low_ess <- parameter_summary[is.finite(parameter_summary$ess_bulk) & parameter_summary$ess_bulk < 100, , drop = FALSE]

diagnostics <- list(
  settings = list(p = p, chains = chains, iter_warmup = iter_warmup,
                  iter_sampling = iter_sampling, thin = thin,
                  retained_draws_per_chain = iter_sampling / thin,
                  adapt_delta = adapt_delta,
                  max_treedepth = max_treedepth, seed = seed),
  data_dimensions = c(t = nrow(y), N = ncol(y), p = p),
  divergences = divergences,
  treedepth_hits = treedepth_hits,
  ebfmi_by_chain = ebfmi,
  parameter_summary = parameter_summary,
  parameters_with_rhat_over_1.01 = bad_rhat,
  parameters_with_bulk_ess_below_100 = low_ess,
  cmdstan_diagnostic_summary = fit$diagnostic_summary()
)

saveRDS(diagnostics, file.path(output_dir, "diagnostics.rds"))
saveRDS(draws, file.path(output_dir, "posterior_draws.rds"))
saveRDS(sampler_draws, file.path(output_dir, "sampler_diagnostics.rds"))
write.csv(parameter_summary, file.path(output_dir, "parameter_summary.csv"), row.names = FALSE)
fit$save_object(file.path(output_dir, "fit.rds"))

cat("\n================ VAR-SV DIAGNOSTICS ================\n")
cat("Divergences: ", divergences, "\n", sep = "")
cat("Max-treedepth hits: ", treedepth_hits, "\n", sep = "")
cat("E-BFMI by chain:\n")
print(ebfmi)
cat("Parameters with R-hat > 1.01: ", nrow(bad_rhat), "\n", sep = "")
cat("Parameters with bulk ESS < 100: ", nrow(low_ess), "\n", sep = "")
cat("Results saved in: ", normalizePath(output_dir), "\n", sep = "")
