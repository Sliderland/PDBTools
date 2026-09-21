# Batch example for PDBEntryBuilder_v2
#
# This example defines two related posterior entries that share one Stan model.
# `run_workflows()` accepts them in one batch, registers the common model once,
# and processes each posterior sequentially. It does not sample in parallel.
#
# The defaults below are a safe preview: sourcing this file constructs the
# specifications and checks the batch interface without writing to PosteriorDB
# or starting Stan sampling. Change the flags deliberately when ready.

source("PDBEntryBuilder_v2.R")

pdb_path <- path.expand("~/Documents/posteriordb")
bayesian_ts_path <- file.path(
    path.expand("~/OneDrive - Uppsala universitet"),
    "Bayesian Time Series",
    "Bayesian Time Series Code"
)
heaps_program_path <- file.path(bayesian_ts_path, "HeapsStanPrograms")

# Execution controls -------------------------------------------------------

register_entries <- TRUE
run_sampling <- TRUE
write_reference_files <- TRUE
overwrite_existing <- TRUE
continue_on_error <- TRUE
save_failed_fits <- FALSE
failed_fit_dir <- path.expand(
    "~/Documents/PDBTools_diagnostics/failed_reference_fits"
)

if (write_reference_files && !run_sampling) {
    stop("`write_reference_files = TRUE` requires `run_sampling = TRUE`.")
}
if (run_sampling && !register_entries) {
    stop("Sampling this example requires `register_entries = TRUE`.")
}

entry <- PDBEntryBuilder$new(pdb_path)

# Shared bibliography entry -----------------------------------------------

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

if (register_entries) {
    # add_bibtex_entry() is duplicate-safe, so this is idempotent by key.
    entry$add_bibtex_entry(heaps_reference)
}

# Construct the related datasets ------------------------------------------

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
        p = 4,
        N = nrow(processed$y),
        y = processed$y,
        es = c(0, 0),
        fs = sqrt(c(0.455, 0.455)),
        gs = c(1.365, 1.365),
        hs = c(0.071175, 0.071175),
        scale_diag = 1,
        scale_offdiag = 0,
        df = dimension + 4,
        grainsize = 25
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

# Define the model once ----------------------------------------------------

statprior_model <- list(
    stan_file = file.path(heaps_program_path, "statpriorPDB.stan"),
    info = list(
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
            "A stationary, exchangeable prior that ensures stationarity in",
            "a VAR model. It maps VAR coefficients to partial",
            "autocorrelations and reverses the mapping."
        ),
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
)

# Shared posterior and sampling settings ----------------------------------

excluded_parameters <- c(
    "topblock",
    "companion",
    "lambdas",
    "lambda_moduli",
    "max_lambda_modulus"
)

# This matches the long-run configuration used for the Heaps models:
# 10 chains x 1,000 retained post-warmup draws = 10,000 reference draws.
sampling_args <- list(
    chains = 10,
    iter = 30000,
    warmup = 10000,
    refresh = 10000,
    thin = 20,
    seed = 123,
    control = list(adapt_delta = 0.95)
)

# Build one entry specification per dataset -------------------------------

entry_specs <- Map(
    function(name, data, info) {
        list(
            data = list(
                data = data,
                info = info
            ),
            posterior = list(
                parameters_exclude = excluded_parameters
                # `dimensions` is intentionally omitted so v2 infers it.
            ),
            reference = list(
                sampling_args = sampling_args,
                comments = paste(
                    "Batch-generated reference draws for",
                    paste(name, "statprior_var", sep = "-")
                )
            )
        )
    },
    names(heaps_data),
    heaps_data,
    heaps_data_info
)
names(entry_specs) <- names(heaps_data)

# Execute the batch --------------------------------------------------------

batch_results <- entry$run_workflows(
    model = statprior_model,
    entries = entry_specs,
    register = register_entries,
    sample = run_sampling,
    write = write_reference_files,
    overwrite = overwrite_existing,
    continue_on_error = continue_on_error,
    save_failed_fits = save_failed_fits,
    failed_fit_dir = failed_fit_dir
)

batch_summary <- entry$summarize_workflow_results(batch_results)
print(batch_summary)

# Individual results remain available for inspection, for example:
batch_results$heaps_small3$status
batch_results$heaps_small3$fit
batch_results$heaps_small3$diagnostics
batch_results$heaps_small3$failed_required_checks
batch_results$heaps_small3$summary_paths

batch_results$heaps_med10$status
batch_results$heaps_med10$fit
batch_results$heaps_med10$diagnostics
batch_results$heaps_med10$failed_required_checks
batch_results$heaps_med10$summary_paths

batch_results$heaps_med20$status
batch_results$heaps_med20$fit
batch_results$heaps_med20$diagnostics
batch_results$heaps_med20$failed_required_checks
batch_results$heaps_med20$summary_paths
