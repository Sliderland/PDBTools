# Batch registration and reference sampling for a canonical VAR/VECM set.
#
# The file is deliberately staged: preflight is the default, and registration,
# sampling, and reference-file writing are separate switches.  Set
# `source_model_root` to the Bayesian Time Series Code directory containing
# the Stan programs, then enable the stages deliberately.

# Configuration ------------------------------------------------------------

source_model_root <- Sys.getenv(
    "BAYESIAN_TS_CODE",
    unset = file.path(
        "/Users/gerpr308/Library/CloudStorage",
        "OneDrive-Uppsalauniversitet",
        "Bayesian Time Series",
        "Bayesian Time Series Code"
    )
)

pdb_path <- Sys.getenv(
    "POSTERIORDB_PATH",
    unset = "/Users/gerpr308/Documents/posteriordb"
)

# Safe default.  This checks paths and prints the proposed Cartesian plan.
preflight_only <- TRUE
register_entries <- FALSE
run_sampling <- FALSE
write_reference_files <- FALSE

overwrite_registration <- FALSE
overwrite_reference_files <- FALSE
skip_completed_references <- TRUE
continue_on_error <- TRUE

# Retain completed but non-passing rstanfit objects outside PosteriorDB for
# later diagnostics. These files can be large; leave disabled unless needed.
save_failed_fits <- FALSE
failed_fit_dir <- file.path(
    path.expand("~"),
    "Documents",
    "PDBTools_diagnostics",
    "failed_reference_fits"
)
save_all_fits <- FALSE
all_fit_dir <- path.expand(
    "~/Documents/PDBTools_diagnostics/all_reference_fits"
)

dataset_sizes_to_run <- c("small3", "med10")
models_to_run <- c(
    "unconstrained_var",
    "minnesota_var",
    "horseshoe_var",
    "vecm_priors",
    "vecm_long_run"
)

# Optional after repairing the current Stan sources.  They remain fully
# described below so they can be enabled without changing the batch logic.
# - VAR_SV currently assigns a Beta prior to a parameter with support [-1, 1].
# - VECM_URCA_test still has flat priors on mu, alpha, and phi and needs
#   data-specific validation before reference sampling.
# models_to_run <- c(models_to_run, "var_sv", "vecm_urca")

var_lag <- 4L
vecm_rank <- 2L

if (write_reference_files && !run_sampling) {
    stop("`write_reference_files = TRUE` requires `run_sampling = TRUE`.")
}
if (preflight_only && any(c(register_entries, run_sampling, write_reference_files))) {
    stop("`preflight_only` cannot be combined with an active workflow stage.")
}

# Model catalogue ----------------------------------------------------------

model_definitions <- list(
    unconstrained_var = list(
        stan_file = file.path(source_model_root, "VARP_as_VAR1_PDB.stan"),
        title = "Unconstrained VAR(p) Model",
        description = paste(
            "A Bayesian vector autoregression with an intercept, unrestricted",
            "autoregressive coefficient matrices, and correlated innovations."
        ),
        keywords = c("VAR", "BVAR", "unconstrained", "correlated innovations"),
        references = "sims1980macroeconomics",
        parameters_exclude = c(
            "topblock", "companion", "lambdas", "lambda_moduli",
            "max_lambda_modulus"
        )
    ),
    minnesota_var = list(
        stan_file = file.path(source_model_root, "MinnesotaVARP.stan"),
        title = "Minnesota-Prior VAR(p) Model",
        description = paste(
            "A Bayesian VAR with a random-walk prior on own first lags and",
            "hierarchical lag and cross-variable shrinkage."
        ),
        keywords = c("VAR", "BVAR", "Minnesota prior", "shrinkage"),
        references = "litterman1986forecasting",
        parameters_exclude = c(
            "companion", "lambdas", "lambda_moduli", "max_lambda_modulus"
        )
    ),
    horseshoe_var = list(
        stan_file = file.path(source_model_root, "HSVAR.stan"),
        title = "Horseshoe-Shrinkage VAR(p) Model",
        description = paste(
            "A Bayesian VAR using global-local horseshoe shrinkage for the",
            "autoregressive coefficients and a correlated innovation covariance."
        ),
        keywords = c("VAR", "BVAR", "horseshoe", "global-local shrinkage"),
        references = "carvalho2010horseshoe",
        parameters_exclude = character()
    ),
    var_sv = list(
        stan_file = file.path(source_model_root, "VAR_SV.stan"),
        title = "VAR(p) with Stochastic Volatility",
        description = paste(
            "A constant-coefficient VAR with time-varying diagonal",
            "innovation standard deviations following autoregressive log-scale",
            "processes."
        ),
        keywords = c("VAR", "BVAR", "stochastic volatility", "heteroskedasticity"),
        references = "primiceri2005time",
        parameters_exclude = character()
    ),
    vecm_priors = list(
        stan_file = file.path(source_model_root, "VECMwPriors.stan"),
        title = "Bayesian Transitory VECM(p)",
        description = paste(
            "A Bayesian vector error-correction model with normalized",
            "cointegration vectors, short-run dynamics, and correlated innovations."
        ),
        keywords = c("VECM", "cointegration", "VAR", "error correction"),
        references = "johansen1995likelihood",
        parameters_exclude = c(
            "A", "topblock", "ect", "short_run", "inLevelsTop",
            "inLevelsCompanion", "lambdas", "lambda_moduli",
            "max_lambda_modulus"
        )
    ),
    vecm_long_run = list(
        stan_file = file.path(source_model_root, "VECMLongRun.stan"),
        title = "Bayesian Long-Run VECM(p)",
        description = paste(
            "A Bayesian VECM whose error-correction term is evaluated at the",
            "long-run lag, with normalized cointegration vectors and correlated",
            "innovations."
        ),
        keywords = c("VECM", "cointegration", "long-run", "VAR"),
        references = "johansen1995likelihood",
        parameters_exclude = c(
            "A", "topblock", "ect", "short_run", "inLevelsTop",
            "inLevelsCompanion", "lambdas", "lambda_moduli",
            "max_lambda_modulus"
        )
    ),
    vecm_urca = list(
        stan_file = file.path(source_model_root, "VECM_URCA_test.stan"),
        title = "Bayesian VECM(p) with Exogenous Deterministic Terms",
        description = paste(
            "A Bayesian VECM with normalized cointegration vectors, short-run",
            "dynamics, correlated innovations, and user-supplied exogenous",
            "deterministic regressors such as seasonal dummies."
        ),
        keywords = c("VECM", "cointegration", "seasonality", "exogenous", "VAR"),
        references = "johansen1995likelihood",
        parameters_exclude = c(
            "A", "topblock", "ect", "short_run", "inLevelsTop",
            "inLevelsCompanion", "lambdas", "lambda_moduli",
            "max_lambda_modulus"
        )
    )
)

unknown_models <- setdiff(models_to_run, names(model_definitions))
if (length(unknown_models)) {
    stop("Unknown model selections: ", paste(unknown_models, collapse = ", "))
}
model_definitions <- model_definitions[models_to_run]

# Bibliography -------------------------------------------------------------

bibliography_entries <- c(
    paste(
    "@article{sims1980macroeconomics,",
    "  title = {Macroeconomics and Reality},",
    "  author = {Sims, Christopher A.},",
    "  journal = {Econometrica}, year = {1980}, volume = {48},",
    "  number = {1}, pages = {1--48}}",
    sep = "\n"
    ),
    paste(
    "@article{litterman1986forecasting,",
    "  title = {Forecasting with {Bayesian} Vector Autoregressions---Five",
    "  Years of Experience}, author = {Litterman, Robert B.},",
    "  Journal = {Journal of Business and Economic Statistics},",
    "  year = {1986}, volume = {4}, number = {1}, pages = {25--38}}",
    sep = "\n"
    ),
    paste(
    "@article{carvalho2010horseshoe,",
    "  title = {The Horseshoe Estimator for Sparse Signals},",
    "  author = {Carvalho, Carlos M. and Polson, Nicholas G. and Scott, James G.},",
    "  journal = {Biometrika}, year = {2010}, volume = {97}, number = {2},",
    "  pages = {465--480}}",
    sep = "\n"
    ),
    paste(
    "@article{primiceri2005time,",
    "  title = {Time Varying Structural Vector Autoregressions and Monetary Policy},",
    "  author = {Primiceri, Giorgio E.}, Journal = {Review of Economic Studies},",
    "  year = {2005}, volume = {72}, number = {3}, pages = {821--852}}",
    sep = "\n"
    ),
    paste(
    "@book{johansen1995likelihood,",
    "  title = {Likelihood-Based Inference in Cointegrated Vector Autoregressive Models},",
    "  author = {Johansen, Soren}, publisher = {Oxford University Press}, year = {1995}}",
    sep = "\n"
    )
)

# Data builders ------------------------------------------------------------

if (!dir.exists(pdb_path)) {
    stop("PDB path does not exist: ", pdb_path)
}
missing_stan <- vapply(
    model_definitions,
    function(x) !file.exists(x$stan_file),
    logical(1)
)
if (any(missing_stan)) {
    stop(
        "Missing Stan source(s): ",
        paste(vapply(model_definitions[missing_stan], `[[`, character(1), "stan_file"), collapse = ", ")
    )
}

dataset_definitions <- list(
    small3 = list(dimension = 3L, size = "Small"),
    med10 = list(dimension = 10L, size = "Medium"),
    med20 = list(dimension = 20L, size = "Large")
)
unknown_sizes <- setdiff(dataset_sizes_to_run, names(dataset_definitions))
if (length(unknown_sizes)) {
    stop("Unknown dataset selections: ", paste(unknown_sizes, collapse = ", "))
}
dataset_definitions <- dataset_definitions[dataset_sizes_to_run]

var_model_names <- c(
    "unconstrained_var", "minnesota_var", "horseshoe_var", "var_sv"
)
selected_var_models <- intersect(var_model_names, names(model_definitions))
heaps_y <- list()
if (length(selected_var_models)) {
    heaps_read <- file.path(source_model_root, "HeapsStanPrograms", "read.R")
    if (!file.exists(heaps_read)) {
        stop("Missing Heaps data reader: ", heaps_read)
    }
    source(heaps_read)
    heaps_y <- lapply(dataset_definitions, function(x) {
        process_data(
            num = x$dimension,
            omit = c(1:2, 197:200),
            Nahead = 40,
            data_dir_path = file.path(source_model_root, "HeapsStanPrograms", "data")
        )$y
    })
}

build_var_data <- function(model_name, y) {
    n <- ncol(y)
    if (model_name == "unconstrained_var") {
        return(list(m = n, p = var_lag, N = nrow(y), y = y))
    }
    if (model_name == "minnesota_var") {
        return(list(
            t = nrow(y), N = n, p = var_lag, d = n + 4,
            Sigma0 = diag(n), lambda0 = 1, theta0 = 1, y = y
        ))
    }
    if (model_name == "horseshoe_var") {
        return(list(t = nrow(y), N = n, p = var_lag, tau_0 = 0.1, y = y))
    }
    if (model_name == "var_sv") {
        return(list(t = nrow(y), N = n, p = var_lag, y = y))
    }
    stop("Not a standard VAR model: ", model_name)
}

make_data_info <- function(name, title, description, keywords, references) {
    list(
        name = name,
        keywords = keywords,
        title = title,
        description = description,
        references = references,
        added_by = "Gerald Press",
        added_date = Sys.Date()
    )
}

data_entries <- list()
for (dataset_key in names(heaps_y)) {
    y <- heaps_y[[dataset_key]]
    dataset_meta <- dataset_definitions[[dataset_key]]
    for (model_name in names(model_definitions)) {
        if (!model_name %in% var_model_names) next
        name <- paste("var_literature", dataset_key, model_name, sep = "_")
        data_entries[[name]] <- list(
            model_name = model_name,
            data = build_var_data(model_name, y),
            info = make_data_info(
                name,
                paste(dataset_meta$size, "Heaps macroeconomic data for", model_name),
                paste(
                    "Quarterly US macroeconomic observations with",
                    dataset_meta$dimension, "variables and model-specific",
                    "Stan data fields."
                ),
                c("macroeconomics", "US", "quarterly", "VAR", model_name),
                model_definitions[[model_name]]$references
            )
        )
    }
}

selected_vecm_models <- intersect(
    c("vecm_priors", "vecm_long_run"),
    names(model_definitions)
)
shared_vecm_data_name <- NULL
if (length(selected_vecm_models)) {
    source("VECMSyntheticData.R")
    vecm_simulation <- generate_vecm_synthetic_data()
    vecm_base <- vecm_simulation$stan_data
    if (vecm_base$p != var_lag || vecm_base$h != vecm_rank) {
        stop("The VECM generator and selected lag/rank settings disagree.")
    }
    shared_vecm_data_name <- "vecm_synthetic5_rank2_p4_v1"
    data_entries[[shared_vecm_data_name]] <- list(
        model_name = selected_vecm_models[[1L]],
        data = vecm_base,
        info = make_data_info(
            shared_vecm_data_name,
            "Five-variable rank-two synthetic cointegrated series",
            paste(
                vecm_base$t,
                "observations from a Gaussian VECM(4) with rank two,",
                "three common stochastic trends, correlated innovations,",
                "and seed 20260917. The first four rows are fixed",
                "conditioning values. Reproduce with VECMSyntheticData.R."
            ),
            c("VECM", "cointegration", "synthetic", "rank 2"),
            "johansen1995likelihood"
        )
    )
    message(
        "Synthetic VECM checks: rank=", vecm_simulation$checks$rank,
        ", unit roots=", vecm_simulation$checks$unit_root_count,
        ", max other root=",
        signif(vecm_simulation$checks$max_other_root_modulus, 4)
    )
}

if ("vecm_urca" %in% names(model_definitions)) {
    vecm_archive_path <- file.path(getwd(), "recovered_vecm_data.rds")
    if (!file.exists(vecm_archive_path)) {
        stop("The optional URCA VECM archive is missing: ", vecm_archive_path)
    }
    denmark_data <- readRDS(vecm_archive_path)$stan_data$vecm_urca_hmc
    name <- "vecm_denmark_seasonal_vecm_urca"
    data_entries[[name]] <- list(
        model_name = "vecm_urca",
        data = denmark_data,
        info = make_data_info(
            name,
            "Denmark data with seasonal regressors for the exogenous VECM",
            paste(
                "Five-variable Denmark data with centered quarterly seasonal",
                "dummies and the data fields required by the exogenous VECM."
            ),
            c("VECM", "cointegration", "Denmark", "seasonality", "exogenous"),
            "johansen1995likelihood"
        )
    )
}

plan <- data.frame(
    data_name = names(data_entries),
    model_name = vapply(data_entries, function(x) x$model_name, character(1)),
    stringsAsFactors = FALSE
)
if (length(selected_vecm_models) > 1L) {
    plan <- rbind(
        plan,
        data.frame(
            data_name = rep(
                shared_vecm_data_name,
                length(selected_vecm_models) - 1L
            ),
            model_name = selected_vecm_models[-1L],
            stringsAsFactors = FALSE
        )
    )
}
rownames(plan) <- NULL
plan$posterior_name <- paste(plan$data_name, plan$model_name, sep = "-")
print(plan)

if (preflight_only) {
    message("Preflight completed. No PDB changes or Stan sampling were run.")
    invisible(plan)
}

# Registration and sampling ------------------------------------------------

if (!preflight_only) {
    source("PDBEntryBuilder_v2.R")
    entry <- PDBEntryBuilder$new(pdb_path, n_threads = 2L)

    if (register_entries) {
        for (bibtex_entry in bibliography_entries) {
            entry$add_bibtex_entry(bibtex_entry)
        }

        for (data_name in names(data_entries)) {
            if (
                data_name %in% entry$search_data(data_name) &&
                    !overwrite_registration
            ) {
                message("Skipping existing data: ", data_name)
                next
            }
            entry$add_data(
                data = data_entries[[data_name]]$data,
                info = data_entries[[data_name]]$info,
                overwrite = overwrite_registration
            )
        }

        model_info <- lapply(names(model_definitions), function(model_name) {
            x <- model_definitions[[model_name]]
            list(
                name = model_name,
                keywords = x$keywords,
                title = x$title,
                description = x$description,
                references = x$references,
                framework = "stan",
                added_by = "Gerald Press",
                added_date = Sys.Date()
            )
        })
        names(model_info) <- names(model_definitions)

        for (model_name in names(model_definitions)) {
            if (
                model_name %in% entry$search_model(model_name) &&
                    !overwrite_registration
            ) {
                message("Skipping existing model: ", model_name)
                next
            }
            entry$add_model_code(
                model_definitions[[model_name]]$stan_file,
                model_info[[model_name]],
                overwrite = overwrite_registration
            )
        }

        for (row in seq_len(nrow(plan))) {
            data_name <- plan$data_name[[row]]
            model_name <- plan$model_name[[row]]
            posterior_name <- plan$posterior_name[[row]]
            if (
                posterior_name %in% entry$search_posterior(posterior_name) &&
                    !overwrite_registration
            ) {
                message("Skipping existing posterior: ", posterior_name)
                next
            }
            posterior <- entry$prepare_posterior(
                data_name = data_name,
                model_name = model_name,
                references = model_definitions[[model_name]]$references,
                parameters_exclude = model_definitions[[model_name]]$parameters_exclude
            )
            entry$add_posterior(
                posterior,
                overwrite = overwrite_registration,
                dry_run = FALSE
            )
        }
    }

    sampling_args <- list(
        chains = 10L,
        iter = 30000L,
        warmup = 10000L,
        thin = 20L,
        refresh = 10000L,
        seed = 123L,
        control = list(adapt_delta = 0.95)
    )
    sampling_specs <- setNames(lapply(seq_len(nrow(plan)), function(i) {
        list(
            posterior_name = plan$posterior_name[[i]],
            posterior = list(
                data_name = plan$data_name[[i]],
                model_name = plan$model_name[[i]]
            ),
            reference = list(
                sampling_args = sampling_args,
                comments = paste(
                    "Reference sampling for the canonical VAR/VECM literature set:",
                    plan$posterior_name[[i]]
                )
            )
        )
    }), plan$posterior_name)

    if (skip_completed_references) {
        present <- vapply(
            names(sampling_specs),
            function(x) x %in% entry$search_reference_draws(x),
            logical(1)
        )
        sampling_specs <- sampling_specs[!present]
    }

    if (run_sampling) {
        if (!length(sampling_specs)) {
            message("No new reference draws selected.")
        } else {
            results <- entry$run_workflows(
                entries = sampling_specs,
                register = FALSE,
                sample = TRUE,
                write = write_reference_files,
                overwrite = overwrite_reference_files,
                continue_on_error = continue_on_error,
                save_failed_fits = save_failed_fits,
                failed_fit_dir = failed_fit_dir,
                save_all_fits = save_all_fits,
                all_fit_dir = all_fit_dir
            )
            print(entry$summarize_workflow_results(results))
        }
    }
}
