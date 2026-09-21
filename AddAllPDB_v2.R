# Add the Heaps VAR posterior grid using PDBEntryBuilder_v2.
#
# The intended grid from addAllPDB.R contains three datasets and four models:
#
#   heaps_small3, heaps_med10, heaps_med20
#       x
#   varp, statprior_var, varp_threaded, statprior_var_threaded
#
# This file therefore describes 12 posteriors. Registration and sampling are
# disabled by default. Enable each stage deliberately after reviewing `plan`.

source("PDBEntryBuilder_v2.R")

pdb_path <- path.expand("~/Documents/posteriordb")
bayesian_ts_path <- file.path(
    path.expand("~/OneDrive - Uppsala universitet"),
    "Bayesian Time Series",
    "Bayesian Time Series Code"
)
heaps_program_path <- file.path(bayesian_ts_path, "HeapsStanPrograms")

# Execution controls -------------------------------------------------------

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
failed_fit_dir <- path.expand(
    "~/Documents/PDBTools_diagnostics/failed_reference_fits"
)

# Use these selectors to stage expensive work. The defaults describe all 12
# posteriors, but no work occurs while the execution flags above are FALSE.
data_to_run <- c("heaps_small3", "heaps_med10", "heaps_med20")
models_to_run <- c(
    "varp",
    "statprior_var",
    "varp_threaded",
    "statprior_var_threaded"
)

if (write_reference_files && !run_sampling) {
    stop("`write_reference_files = TRUE` requires `run_sampling = TRUE`.")
}

entry <- PDBEntryBuilder$new(
    path = pdb_path,
    n_threads = 2L
)
rstan::rstan_options(threads_per_chain = 2L)

# Shared bibliography ------------------------------------------------------

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
    "  url = {https://doi.org/10.1080/10618600.2022.2079648}",
    "}",
    sep = "\n"
)

# Data definitions ---------------------------------------------------------

source(file.path(heaps_program_path, "read.R"))

build_heaps_data <- function(dimension) {
    processed <- process_data(
        num = dimension,
        omit = c(1:2, 1970:200),
        Nahead = 40,
        data_dir_path = file.path(heaps_program_path, "data")
    )

    list(
        m = dimension,
        p = 4L,
        N = nrow(processed$y),
        y = processed$y,
        es = c(0, 0),
        fs = sqrt(c(0.455, 0.455)),
        gs = c(1.365, 1.365),
        hs = c(0.071175, 0.071175),
        scale_diag = 1,
        scale_offdiag = 0,
        df = dimension + 4,
        grainsize = 25L
    )
}

make_heaps_data_info <- function(name, dimension, size) {
    list(
        name = name,
        keywords = c(
            "United States",
            "US",
            "Economics",
            "Macroeconomics",
            "Vector Autoregression",
            "BVAR",
            "VAR",
            paste(size, "VAR"),
            size,
            "Quarterly"
        ),
        title = paste0(
            size,
            " Vector Autoregressive Data (n = ",
            dimension,
            ", p = 4)"
        ),
        description = paste0(
            "Quarterly US macroeconomic data. A ",
            tolower(size),
            "-sized, pre-transformed vector autoregressive dataset with ",
            dimension,
            " variables and 4 lags."
        ),
        urls = paste0(
            "https://tandf.figshare.com/articles/dataset/",
            "Enforcing_stationarity_through_the_prior_in_",
            "vector_autoregressions/19831250/2"
        ),
        references = "heaps2023stationary",
        added_date = Sys.Date(),
        added_by = "Gerald Press"
    )
}

data_definitions <- list(
    heaps_small3 = list(dimension = 3L, size = "Small"),
    heaps_med10 = list(dimension = 10L, size = "Medium"),
    heaps_med20 = list(dimension = 20L, size = "Large")
)

# Model definitions --------------------------------------------------------

make_model_info <- function(name, stationary, threaded) {
    model_type <- if (stationary) {
        "Stationary Exchangeable"
    } else {
        "Unconstrained"
    }
    threading <- if (threaded) " Threaded" else ""

    list(
        name = name,
        keywords = c(
            if (stationary) c("exchangeable", "stationary") else {
                "unconstrained"
            },
            if (threaded) "threaded" else NULL,
            "VAR",
            "Vector Autoregressive",
            "Vector Autoregression",
            "BVAR"
        ),
        title = paste0(model_type, threading, " VAR(p) Model"),
        description = if (stationary) {
            paste(
                "A stationary, exchangeable prior that ensures stationarity",
                "in a VAR model. It maps VAR coefficients to partial",
                "autocorrelations and reverses the mapping."
            )
        } else {
            "An unconstrained VAR(p) model."
        },
        urls = paste0(
            "https://tandf.figshare.com/articles/dataset/",
            "Enforcing_stationarity_through_the_prior_in_",
            "vector_autoregressions/19831250/2"
        ),
        references = "heaps2023stationary",
        framework = "stan",
        added_by = "Gerald Press",
        added_date = Sys.Date()
    )
}

generated_quantity_exclusions <- c(
    "topblock",
    "companion",
    "lambdas",
    "lambda_moduli",
    "max_lambda_modulus"
)

model_definitions <- list(
    varp = list(
        stan_file = file.path(bayesian_ts_path, "VARP_as_VAR1_PDB.stan"),
        stationary = FALSE,
        threaded = FALSE,
        parameters_exclude = generated_quantity_exclusions
    ),
    statprior_var = list(
        stan_file = file.path(heaps_program_path, "statpriorPDB.stan"),
        stationary = TRUE,
        threaded = FALSE,
        parameters_exclude = generated_quantity_exclusions
    ),
    varp_threaded = list(
        stan_file = file.path(bayesian_ts_path, "var_threaded_grain.stan"),
        stationary = FALSE,
        threaded = TRUE,
        parameters_exclude = generated_quantity_exclusions
    ),
    statprior_var_threaded = list(
        stan_file = file.path(bayesian_ts_path, "statprior_threaded.stan"),
        stationary = TRUE,
        threaded = TRUE,
        parameters_exclude = generated_quantity_exclusions
    )
)

unknown_data <- setdiff(data_to_run, names(data_definitions))
unknown_models <- setdiff(models_to_run, names(model_definitions))
if (length(unknown_data) > 0L) {
    stop("Unknown data selections: ", paste(unknown_data, collapse = ", "))
}
if (length(unknown_models) > 0L) {
    stop("Unknown model selections: ", paste(unknown_models, collapse = ", "))
}

data_definitions <- data_definitions[data_to_run]
model_definitions <- model_definitions[models_to_run]

heaps_data <- lapply(data_definitions, function(definition) {
    build_heaps_data(definition$dimension)
})
heaps_data_info <- Map(
    function(name, definition) {
        make_heaps_data_info(
            name = name,
            dimension = definition$dimension,
            size = definition$size
        )
    },
    names(data_definitions),
    data_definitions
)

model_info <- Map(
    function(name, definition) {
        make_model_info(
            name = name,
            stationary = definition$stationary,
            threaded = definition$threaded
        )
    },
    names(model_definitions),
    model_definitions
)

plan <- expand.grid(
    data_name = names(data_definitions),
    model_name = names(model_definitions),
    stringsAsFactors = FALSE
)
plan$posterior_name <- paste(
    plan$data_name,
    plan$model_name,
    sep = "-"
)
print(plan)

# Sampling configuration --------------------------------------------------

# 10 chains x 1,000 retained post-warmup draws = 10,000 reference draws.
sampling_args <- list(
    chains = 10L,
    iter = 30000L,
    warmup = 10000L,
    refresh = 10000L,
    thin = 20L,
    seed = 123L,
    control = list(adapt_delta = 0.95)
)

registration_results <- list(
    bibliography = NULL,
    data = list(),
    models = list(),
    posteriors = list()
)

attempt <- function(label, expression) {
    tryCatch(
        list(status = "completed", value = force(expression), error = NULL),
        error = function(error) {
            message(label, " failed: ", conditionMessage(error))
            list(
                status = "error",
                value = NULL,
                error = conditionMessage(error)
            )
        }
    )
}

exists_exactly <- function(matches, name) {
    name %in% matches
}

# Register shared objects and posterior definitions -----------------------

if (register_entries) {
    registration_results$bibliography <- attempt(
        "Bibliography registration",
        entry$add_bibtex_entry(heaps_reference)
    )

    for (data_name in names(data_definitions)) {
        already_exists <- exists_exactly(
            entry$search_data(data_name),
            data_name
        )
        registration_results$data[[data_name]] <- if (
            already_exists && !overwrite_registration
        ) {
            list(status = "skipped_existing", value = NULL, error = NULL)
        } else {
            attempt(
                paste("Data registration", data_name),
                entry$add_data(
                    data = heaps_data[[data_name]],
                    info = heaps_data_info[[data_name]],
                    overwrite = overwrite_registration
                )
            )
        }
    }

    for (model_name in names(model_definitions)) {
        definition <- model_definitions[[model_name]]
        already_exists <- exists_exactly(
            entry$search_model(model_name),
            model_name
        )
        registration_results$models[[model_name]] <- if (
            already_exists && !overwrite_registration
        ) {
            list(status = "skipped_existing", value = NULL, error = NULL)
        } else {
            attempt(
                paste("Model registration", model_name),
                entry$add_model_code(
                    stan_file = definition$stan_file,
                    info = model_info[[model_name]],
                    overwrite = overwrite_registration
                )
            )
        }
    }

    for (row in seq_len(nrow(plan))) {
        data_name <- plan$data_name[[row]]
        model_name <- plan$model_name[[row]]
        posterior_name <- plan$posterior_name[[row]]
        already_exists <- exists_exactly(
            entry$search_posterior(posterior_name),
            posterior_name
        )

        registration_results$posteriors[[posterior_name]] <- if (
            already_exists && !overwrite_registration
        ) {
            list(status = "skipped_existing", value = NULL, error = NULL)
        } else {
            attempt(
                paste("Posterior registration", posterior_name),
                {
                    posterior <- entry$prepare_posterior(
                        data_name = data_name,
                        model_name = model_name,
                        references = "heaps2023stationary",
                        parameters_exclude = model_definitions[[model_name]]$parameters_exclude
                    )
                    entry$add_posterior(
                        spec = posterior,
                        overwrite = overwrite_registration,
                        dry_run = FALSE
                    )
                }
            )
        }
    }
}

# Build the sampling batch -------------------------------------------------

sampling_specs <- setNames(
    lapply(seq_len(nrow(plan)), function(row) {
        data_name <- plan$data_name[[row]]
        model_name <- plan$model_name[[row]]
        posterior_name <- plan$posterior_name[[row]]

        list(
            posterior_name = posterior_name,
            posterior = list(
                data_name = data_name,
                model_name = model_name
            ),
            reference = list(
                sampling_args = sampling_args,
                comments = paste(
                    "Sequential batch reference draws for",
                    posterior_name
                )
            )
        )
    }),
    plan$posterior_name
)

skipped_reference_names <- character()
if (skip_completed_references) {
    completed <- vapply(
        names(sampling_specs),
        function(posterior_name) {
            exists_exactly(
                entry$search_reference_draws(posterior_name),
                posterior_name
            )
        },
        logical(1)
    )
    skipped_reference_names <- names(sampling_specs)[completed]
    sampling_specs <- sampling_specs[!completed]
}

if (length(skipped_reference_names) > 0L) {
    message(
        "Skipping existing reference draws: ",
        paste(skipped_reference_names, collapse = ", ")
    )
}

batch_results <- NULL
batch_summary <- NULL

if (run_sampling) {
    if (length(sampling_specs) == 0L) {
        message("Every selected reference posterior already exists.")
    } else {
        batch_results <- entry$run_workflows(
            model = NULL,
            entries = sampling_specs,
            register = FALSE,
            sample = TRUE,
            write = write_reference_files,
            overwrite = overwrite_reference_files,
            continue_on_error = continue_on_error,
            save_failed_fits = save_failed_fits,
            failed_fit_dir = failed_fit_dir
        )
        batch_summary <- entry$summarize_workflow_results(batch_results)
        print(batch_summary)
    }
} else {
    message(
        "Sampling is disabled. Set `run_sampling <- TRUE` after registering ",
        "and reviewing the plan."
    )
}

# Useful inspection commands:
# registration_results
# batch_summary
# batch_results[["heaps_small3-statprior_var"]]$status
# batch_results[["heaps_small3-statprior_var"]]$error
# batch_results[["heaps_small3-statprior_var"]]$failed_required_checks
