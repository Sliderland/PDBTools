# Integration workflow for heaps3-statprior_var
#
# This reproduces the data construction, Stan program, excluded parameters,
# and sampling arguments used in addStatPriorVAR.R. Sampling is deliberately
# kept separate from registration so either stage can be rerun on its own.

source("PDBEntryBuilder_v1.R")

pdb_path <- "/Users/gerpr308/Documents/posteriordb"
bayesian_ts_path <- paste0(
    "/Users/gerpr308/OneDrive - Uppsala universitet/",
    "Bayesian Time Series/Bayesian Time Series Code"
)

register_entries <- TRUE
run_sampling <- TRUE
write_reference_files <- TRUE
overwrite_test_entries <- TRUE
n_ess_failures_to_show <- 10L

if (write_reference_files && !run_sampling) {
    stop("`write_reference_files = TRUE` requires `run_sampling = TRUE`.")
}

entry <- PDBEntryBuilder$new(pdb_path)

heaps_reference_key <- "heaps2023stationary"
heaps_reference <- paste(
    "@article{heaps2023stationary,",
    "  title = {Enforcing {Stationarity} through the {Prior} in {Vector} {Autoregressions}},",
    "  author = {Heaps, Sarah E.},",
    "  journal = {Journal of Computational and Graphical Statistics},",
    "  year = {2023},",
    "  volume = {32},",
    "  number = {1},",
    "  pages = {74--83},",
    "  doi = {10.1080/10618600.2022.2079648},",
    "  url = {https://doi.org/10.1080/10618600.2022.2079648},",
    "  issn = {1061-8600},",
    "  publisher = {Taylor \\& Francis},",
    "  keywords = {Partial autocorrelation matrix, Stan, Unconstrained reparameterization, Vector autoregressive model}",
    "}",
    sep = "\n"
)

references_path <- file.path(
    pdb_path,
    "posterior_database",
    "bibliography",
    "references.bib"
)

reference_is_registered <- file.exists(references_path) &&
    any(grepl(
        paste0("{", heaps_reference_key, ","),
        readLines(references_path, warn = FALSE),
        fixed = TRUE
    ))

if (register_entries && !reference_is_registered) {
    message("Adding bibliography entry: ", heaps_reference_key)
    entry$add_bibtex_entry(heaps_reference)
}

# process_data() and its Stock-Watson preprocessing helpers are defined here.
source(file.path(bayesian_ts_path, "HeapsStanPrograms", "read.R"))

data_builder_heaps <- function(my_m, my_path) {
    my_omit <- c(1:2, 1970:200)
    my_Nahead <- 40
    yraw_trim <- process_data(
        my_m,
        my_omit,
        my_Nahead,
        file.path(my_path, "data")
    )
    yraw_trim_hb <- yraw_trim$y

    my_p <- 4
    list(
        m = my_m,
        p = my_p,
        N = nrow(yraw_trim_hb),
        y = yraw_trim_hb,
        es = c(0, 0),
        fs = sqrt(c(0.455, 0.455)),
        gs = c(1.365, 1.365),
        hs = c(0.071175, 0.071175),
        scale_diag = 1,
        scale_offdiag = 0,
        df = my_m + 4,
        grainsize = 25
    )
}

heaps_program_path <- file.path(bayesian_ts_path, "HeapsStanPrograms")
heaps3_data <- data_builder_heaps(3, heaps_program_path)

heaps3_info <- list(
    name = "heaps3",
    keywords = c(
        "United States",
        "US",
        "Economics",
        "Macroeconomy",
        "Macroeconomics",
        "Vector Autoregression",
        "BVAR",
        "VAR",
        "Small VAR",
        "Small",
        "Quarterly"
    ),
    title = "Small Vector Autoregressive Data (n = 3, p = 4)",
    description = paste0(
        "Quarterly US Macro-economic data. A small, pre-transformed, ",
        "vector autoregressive dataset for testing Vector Autoregressive ",
        "Models. 3 Variables, 4 Lags."
    ),
    urls = paste0(
        "https://tandf.figshare.com/articles/dataset/",
        "Enforcing_stationarity_through_the_prior_in_vector_autoregressions/",
        "19831250/2"
    ),
    references = "heaps2023stationary",
    added_date = Sys.Date(),
    added_by = "Gerald Press"
)

statprior_var_info <- list(
    name = "statprior_var",
    keywords = c(
        "exchangeable",
        "stationary",
        "VAR",
        "Vector Autoregressive",
        "Vector Autoregression",
        "BVAR"
    ),
    title = "Stationary Exchangeable VAR(p) Model",
    description = paste(
        "A stationary, exchangeable prior that ensures stationarity in a",
        "VAR model. Maps VAR coefficients to partial autocorrelations and",
        "reverses the mapping."
    ),
    urls = paste0(
        "https://tandf.figshare.com/articles/dataset/",
        "Enforcing_stationarity_through_the_prior_in_vector_autoregressions/",
        "19831250/2"
    ),
    references = "heaps2023stationary",
    framework = "stan",
    added_by = "Gerald Press",
    added_date = Sys.Date()
)

posterior_name <- "heaps3-statprior_var"
excluded_params <- c(
    "topblock",
    "companion",
    "lambdas",
    "lambda_moduli",
    "max_lambda_modulus"
)

if (register_entries) {
    entry$add_data(
        data = heaps3_data,
        info = heaps3_info,
        overwrite = overwrite_test_entries
    )

    entry$add_model_code(
        stan_file = file.path(heaps_program_path, "statpriorPDB.stan"),
        info = statprior_var_info,
        overwrite = overwrite_test_entries
    )

    posterior <- entry$prepare_posterior(
        data_name = "heaps3",
        model_name = "statprior_var",
        parameters_exclude = excluded_params
    )

    entry$add_posterior(
        spec = posterior,
        overwrite = overwrite_test_entries,
        dry_run = FALSE
    )
}

# Exact method arguments used for heaps3-statprior_var in addStatPriorVAR.R:
# 10 chains x 1,000 retained post-warmup draws = 10,000 total draws.
sampling_args <- list(
    chains = 10,
    iter = 30000,
    warmup = 10000,
    refresh = 10000,
    thin = 20,
    seed = 123,
    control = list(adapt_delta = 0.9)
)

# heaps3_statprior_var <- rstan::stan(
#     file = file.path(heaps_program_path, "statpriorPDB.stan"),
#     data = heaps3_data,
#     chains = 10,
#     iter = 30000,
#     warmup = 10000,
#     refresh = 10000,
#     thin = 20,
#     seed = 123,
#     control = list(adapt_delta = 0.9)
# )
# rstan::check_divergences(heaps3_statprior_var)
# diag_summ <- rstan::get_sampler_params(heaps3_statprior_var, inc_warmup = FALSE)

reference_result <- NULL

if (run_sampling) {
    message("Starting sampling for ", posterior_name, ".")

    fit <- entry$compute_reference_draws(
        posterior_name = posterior_name,
        sampling_args = sampling_args,
        comments = paste(
            "Reference draws generated using the configuration from",
            "addStatPriorVAR.R."
        ),
        auto_check = FALSE,
        write = FALSE,
        overwrite = overwrite_test_entries
    )

    fit_info <- entry$get_reference_info(fit)
    divergences_by_chain <- fit_info$diagnostics$divergent_transitions
    names(divergences_by_chain) <- paste0(
        "chain_",
        seq_along(divergences_by_chain)
    )
    total_divergences <- sum(divergences_by_chain)

    message(
        posterior_name,
        " finished with ",
        total_divergences,
        " divergent transitions."
    )
    message(
        "Divergences by chain: ",
        paste(
            names(divergences_by_chain),
            divergences_by_chain,
            sep = "=",
            collapse = ", "
        )
    )

    checks <- entry$get_checks_from_stanfit(fit)
    failed_checks <- names(checks)[!vapply(checks, isTRUE, logical(1))]
    required_checks <- c(
        "ndraws_is_10k",
        "nchains_is_gte_4",
        "r_hat_below_1_01",
        "efmi_above_0_2",
        "abs_mean_lag1_ac_below_0_05"
    )
    failed_required_checks <- intersect(failed_checks, required_checks)

    if (!isTRUE(checks$ess_within_bounds)) {
        ess_failures <- entry$get_ess_bounds_failures(fit_info$diagnostics)
        shown_names <- head(ess_failures$any, n_ess_failures_to_show)
        message(
            "Legacy ESS bounds check did not pass: ",
            ess_failures$total_count,
            " unique variables were outside the bounds (bulk: ",
            ess_failures$bulk_count,
            ", tail: ",
            ess_failures$tail_count,
            ")."
        )
        message(
            "First ",
            length(shown_names),
            " affected variables: ",
            paste(shown_names, collapse = ", ")
        )
        message(
            "This deprecated check is recorded but does not prevent writing."
        )
    }
    info_path <- NULL
    draws_path <- NULL
    summary_paths <- NULL

    if (length(failed_required_checks) > 0L) {
        message(
            posterior_name,
            " was not written because these required checks failed: ",
            paste(failed_required_checks, collapse = ", ")
        )
    } else {
        fit <- entry$check_draws_from_stanfit(fit)

        if (write_reference_files) {
            info_path <- entry$write_rpi_from_stan_fit(
                fit,
                overwrite = overwrite_test_entries,
                verify = TRUE
            )
            draws_path <- entry$write_rpd_from_stan_fit(
                fit,
                overwrite = overwrite_test_entries,
                verify = TRUE
            )
            summary_paths <- entry$write_summary_statistics_from_stan_fit(
                fit,
                overwrite = overwrite_test_entries,
                verify = TRUE
            )
            entry$verify_reference_files(fit)
        }
    }

    reference_result <- list(
        posterior_name = posterior_name,
        fit = fit,
        diagnostics = fit_info$diagnostics,
        divergences_by_chain = divergences_by_chain,
        total_divergences = total_divergences,
        checks = checks,
        failed_checks = failed_checks,
        failed_required_checks = failed_required_checks,
        info_path = info_path,
        draws_path = draws_path,
        summary_paths = summary_paths
    )
}
