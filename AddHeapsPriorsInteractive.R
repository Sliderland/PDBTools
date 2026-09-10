# Add Heaps VAR posteriors in computationally staged groups.
#
# Work classes:
# - easy: 3- and 10-variable VARs, excluding statinvert_varma
# - hard: every statinvert_varma model and every 20-variable model
#
# `workflow_mode` controls which group is considered:
# - "easy": run only easy jobs; never prompt for hard jobs
# - "hard": run only hard jobs, after confirmation
# - "all": run easy jobs first, then ask before hard jobs
#
# Declining the hard-job prompt is a successful exit. Results from completed
# easy jobs remain written to PosteriorDB and available in `easy_results`.

source("PDBEntryBuilder_v2.R")
source(file.path("HeapsStanPrograms", "read.R"))

# User controls -----------------------------------------------------------

workflow_mode <- "easy"

# NULL prompts interactively. Set TRUE to approve hard jobs without a prompt,
# or FALSE to decline them without a prompt (useful for non-interactive runs).
hard_confirmation <- NULL

register_entries <- TRUE
write_reference_files <- TRUE
overwrite_registration <- FALSE
overwrite_reference_files <- FALSE
skip_completed_references <- TRUE
continue_on_error <- TRUE

pdb_path <- path.expand("~/Documents/posteriordb")
heaps_program_path <- normalizePath("HeapsStanPrograms", mustWork = TRUE)

# Four chains and no thinning retain exactly 10,000 post-warmup draws while
# computing 14,000 rather than 300,000 total HMC transitions.
easy_sampling_args <- list(
    chains = 4L,
    iter = 3500L,
    warmup = 1000L,
    thin = 1L,
    refresh = 500L,
    seed = 123L,
    control = list(adapt_delta = 0.9)
)

# Hard jobs use the same efficient baseline. Increase adapt_delta for a
# particular retry only when its diagnostics show divergences.
hard_sampling_args <- easy_sampling_args

workflow_mode <- match.arg(workflow_mode, c("easy", "hard", "all"))
entry <- PDBEntryBuilder$new(pdb_path)

# Shared specifications ---------------------------------------------------

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

dataset_definitions <- list(
    small3 = list(dimension = 3L, size = "Small"),
    med10 = list(dimension = 10L, size = "Medium"),
    med20 = list(dimension = 20L, size = "Large")
)

model_definitions <- list(
    statprior_var = list(
        stan_file = file.path(heaps_program_path, "statpriorPDB.stan"),
        title = "Stationary Exchangeable VAR(p) Model",
        description = paste(
            "A stationary VAR using the exchangeable partial-autocorrelation",
            "prior described in Section 3.2 of Heaps (2023)."
        ),
        keywords = c("stationary", "exchangeable", "VAR", "BVAR")
    ),
    semiconj_var = list(
        stan_file = file.path(heaps_program_path, "semiconj.stan"),
        title = "Semi-Conjugate VAR(p) Model",
        description = paste(
            "A zero-mean VAR with semi-conjugate priors for the",
            "autoregressive coefficients and innovation covariance."
        ),
        keywords = c("semi-conjugate", "VAR", "BVAR", "comparison")
    ),
    statrml_var = list(
        stan_file = file.path(heaps_program_path, "statrmlprior.stan"),
        title = "Stationary Roy-Reparameterized VAR(p) Model",
        description = paste(
            "A stationary VAR with the vague prior based on the Roy et al.",
            "reparameterization described by Heaps (2023)."
        ),
        keywords = c("stationary", "Roy", "reparameterization", "VAR")
    ),
    statinvert_varma = list(
        stan_file = file.path(heaps_program_path, "statinvertprior.stan"),
        title = "Stationary and Invertible VARMA(p, q) Model",
        description = paste(
            "A stationary and invertible VARMA model with exchangeable",
            "priors for the AR and MA partial-autocorrelation parameters."
        ),
        keywords = c("stationary", "invertible", "VARMA", "exchangeable")
    ),
    ansleykohn_var = list(
        stan_file = file.path(heaps_program_path, "ansleykohn.stan"),
        title = "Stationary Ansley-Kohn VAR(p) Model",
        description = paste(
            "A stationary VAR using the Ansley and Kohn (1986)",
            "partial-autocorrelation reparameterization and uniform priors."
        ),
        keywords = c("stationary", "Ansley-Kohn", "VAR", "uniform prior")
    )
)

prior_defaults <- list(
    p = 4L,
    q = 2L,
    es = c(0, 0),
    fs = sqrt(c(0.455, 0.455)),
    gs = c(1.365, 1.365),
    hs = c(0.071175, 0.071175),
    scale_diag = 1,
    scale_offdiag = 0,
    m_diag = rep(1, 4L),
    s_diag = rep(10, 4L),
    m_offdiag = rep(0, 4L),
    s_offdiag = rep(1, 4L)
)

expected_data_fields <- list(
    statprior_var = c(
        "m", "p", "N", "y", "es", "fs", "gs", "hs",
        "scale_diag", "scale_offdiag", "df"
    ),
    semiconj_var = c(
        "m", "p", "N", "y", "m_diag", "s_diag", "m_offdiag",
        "s_offdiag", "scale_diag", "scale_offdiag", "df"
    ),
    statrml_var = c("m", "p", "N", "y", "df"),
    statinvert_varma = c(
        "m", "p", "q", "N", "y", "es", "fs", "gs", "hs",
        "scale_diag", "scale_offdiag", "df"
    ),
    ansleykohn_var = c("m", "p", "N", "y")
)

build_observations <- function(dimension) {
    process_data(
        num = dimension,
        omit = c(1:2, 197:200),
        Nahead = 40,
        data_dir_path = file.path(heaps_program_path, "data")
    )$y
}

build_model_data <- function(model_name, y) {
    base <- list(
        m = ncol(y),
        p = prior_defaults$p,
        N = nrow(y),
        y = y
    )
    covariance_prior <- list(
        scale_diag = prior_defaults$scale_diag,
        scale_offdiag = prior_defaults$scale_offdiag,
        df = ncol(y) + 4
    )

    data <- switch(
        model_name,
        statprior_var = c(
            base,
            prior_defaults[c("es", "fs", "gs", "hs")],
            covariance_prior
        ),
        semiconj_var = c(
            base,
            prior_defaults[c("m_diag", "s_diag", "m_offdiag", "s_offdiag")],
            covariance_prior
        ),
        statrml_var = c(base, list(df = ncol(y) + 4)),
        statinvert_varma = c(
            base,
            list(
                q = prior_defaults$q,
                es = rbind(prior_defaults$es, prior_defaults$es),
                fs = rbind(prior_defaults$fs, prior_defaults$fs),
                gs = rbind(prior_defaults$gs, prior_defaults$gs),
                hs = rbind(prior_defaults$hs, prior_defaults$hs)
            ),
            covariance_prior
        ),
        ansleykohn_var = base,
        stop("Unknown model: ", model_name, call. = FALSE)
    )

    expected <- expected_data_fields[[model_name]]
    missing_fields <- setdiff(expected, names(data))
    unnecessary_fields <- setdiff(names(data), expected)
    if (length(missing_fields) > 0L || length(unnecessary_fields) > 0L) {
        stop(
            "Invalid model-specific data fields for `", model_name, "`.",
            call. = FALSE
        )
    }
    data
}

generated_quantity_exclusions <- c(
    "topblock",
    "companion",
    "lambdas",
    "lambda_moduli",
    "max_lambda_modulus"
)

# Build the complete plan without registering or sampling -----------------

observations <- lapply(dataset_definitions, function(definition) {
    build_observations(definition$dimension)
})

workflow_entries <- list()
for (dataset_key in names(dataset_definitions)) {
    dataset_definition <- dataset_definitions[[dataset_key]]
    y <- observations[[dataset_key]]

    for (model_name in names(model_definitions)) {
        data_name <- if (identical(model_name, "statprior_var")) {
            paste0("heaps_", dataset_key)
        } else {
            paste0("heaps_", dataset_key, "_", model_name)
        }
        posterior_name <- paste(data_name, model_name, sep = "-")
        is_hard <- identical(model_name, "statinvert_varma") ||
            identical(dataset_key, "med20")

        workflow_entries[[posterior_name]] <- list(
            class = if (is_hard) "hard" else "easy",
            dataset_key = dataset_key,
            model_name = model_name,
            posterior_name = posterior_name,
            data = list(
                data = build_model_data(model_name, y),
                info = list(
                    name = data_name,
                    keywords = c(
                        "United States", "Macroeconomics",
                        "Vector Autoregression", dataset_definition$size,
                        model_name
                    ),
                    title = paste(
                        dataset_definition$size,
                        "Heaps macroeconomic data for",
                        model_name
                    ),
                    description = paste0(
                        "Quarterly US macroeconomic observations with ",
                        dataset_definition$dimension,
                        " variables and model-specific hyperparameters for ",
                        model_name,
                        "."
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
            ),
            posterior = list(
                data_name = data_name,
                model_name = model_name,
                references = "heaps2023stationary",
                parameters_exclude = generated_quantity_exclusions
            )
        )
    }
}

workflow_plan <- data.frame(
    posterior_name = names(workflow_entries),
    class = vapply(workflow_entries, `[[`, character(1), "class"),
    dataset = vapply(workflow_entries, `[[`, character(1), "dataset_key"),
    model = vapply(workflow_entries, `[[`, character(1), "model_name"),
    row.names = NULL
)
print(workflow_plan)

# Registration and execution helpers -------------------------------------

attempt <- function(label, expression) {
    tryCatch(
        list(status = "completed", value = force(expression), error = NULL),
        error = function(error) {
            message(label, " failed: ", conditionMessage(error))
            list(status = "error", value = NULL, error = conditionMessage(error))
        }
    )
}

register_workflow_entries <- function(entries) {
    results <- list(bibliography = NULL, data = list(), models = list(),
                    posteriors = list())
    if (!register_entries) {
        return(results)
    }

    results$bibliography <- attempt(
        "Bibliography registration",
        entry$add_bibtex_entry(heaps_reference)
    )

    unique_models <- unique(vapply(entries, `[[`, character(1), "model_name"))
    for (model_name in unique_models) {
        definition <- model_definitions[[model_name]]
        results$models[[model_name]] <- attempt(
            paste("Model registration", model_name),
            entry$add_model_code(
                stan_file = definition$stan_file,
                info = list(
                    name = model_name,
                    keywords = definition$keywords,
                    title = definition$title,
                    description = definition$description,
                    urls = "https://doi.org/10.1080/10618600.2022.2079648",
                    references = "heaps2023stationary",
                    framework = "stan",
                    added_by = "Gerald Press",
                    added_date = Sys.Date()
                ),
                overwrite = overwrite_registration
            )
        )
    }

    for (name in names(entries)) {
        spec <- entries[[name]]
        data_name <- spec$data$info$name
        results$data[[data_name]] <- attempt(
            paste("Data registration", data_name),
            entry$add_data(
                data = spec$data$data,
                info = spec$data$info,
                overwrite = overwrite_registration
            )
        )
        results$posteriors[[name]] <- attempt(
            paste("Posterior registration", name),
            {
                posterior <- do.call(entry$prepare_posterior, spec$posterior)
                entry$add_posterior(
                    spec = posterior,
                    overwrite = overwrite_registration,
                    dry_run = FALSE
                )
            }
        )
    }
    results
}

run_workflow_group <- function(entries, sampling_args, label) {
    if (length(entries) == 0L) {
        message("No ", label, " workflows were selected.")
        return(list(registration = NULL, sampling = NULL, summary = NULL))
    }

    message("Registering ", length(entries), " ", label, " workflows.")
    registration <- register_workflow_entries(entries)

    sampling_entries <- lapply(entries, function(spec) {
        list(
            posterior_name = spec$posterior_name,
            posterior = spec$posterior[c("data_name", "model_name")],
            reference = list(
                sampling_args = sampling_args,
                comments = paste("Staged reference sampling for", spec$posterior_name)
            )
        )
    })

    if (skip_completed_references) {
        completed <- vapply(names(sampling_entries), function(name) {
            name %in% entry$search_reference_draws(name)
        }, logical(1))
        if (any(completed)) {
            message(
                "Skipping existing reference draws: ",
                paste(names(sampling_entries)[completed], collapse = ", ")
            )
            sampling_entries <- sampling_entries[!completed]
        }
    }

    if (length(sampling_entries) == 0L) {
        return(list(registration = registration, sampling = list(),
                    summary = data.frame()))
    }

    sampling <- entry$run_workflows(
        model = NULL,
        entries = sampling_entries,
        register = FALSE,
        sample = TRUE,
        write = write_reference_files,
        overwrite = overwrite_reference_files,
        continue_on_error = continue_on_error
    )
    summary <- entry$summarize_workflow_results(sampling)
    print(summary)
    list(registration = registration, sampling = sampling, summary = summary)
}

confirm_hard_jobs <- function(number_of_jobs) {
    if (!is.null(hard_confirmation)) {
        if (!is.logical(hard_confirmation) || length(hard_confirmation) != 1L ||
            is.na(hard_confirmation)) {
            stop("`hard_confirmation` must be NULL, TRUE, or FALSE.", call. = FALSE)
        }
        return(hard_confirmation)
    }
    if (!interactive()) {
        message(
            "Hard sampling was not approved: this is a non-interactive session ",
            "and `hard_confirmation` is NULL."
        )
        return(FALSE)
    }

    answer <- readline(paste0(
        "Run ", number_of_jobs,
        " computationally expensive workflow(s) now? [y/N]: "
    ))
    tolower(trimws(answer)) %in% c("y", "yes")
}

# Execute selected stages -------------------------------------------------

easy_entries <- workflow_entries[
    vapply(workflow_entries, `[[`, character(1), "class") == "easy"
]
hard_entries <- workflow_entries[
    vapply(workflow_entries, `[[`, character(1), "class") == "hard"
]

easy_results <- NULL
hard_results <- NULL
hard_sampling_approved <- FALSE

if (workflow_mode %in% c("easy", "all")) {
    easy_results <- run_workflow_group(
        easy_entries,
        sampling_args = easy_sampling_args,
        label = "easy"
    )
}

if (workflow_mode %in% c("hard", "all")) {
    message("Hard workflows:")
    print(workflow_plan[workflow_plan$class == "hard", , drop = FALSE])
    hard_sampling_approved <- confirm_hard_jobs(length(hard_entries))
    if (hard_sampling_approved) {
        hard_results <- run_workflow_group(
            hard_entries,
            sampling_args = hard_sampling_args,
            label = "hard"
        )
    } else {
        message(
            "Hard sampling declined. All completed easy workflows remain saved; ",
            "the script finished successfully."
        )
    }
}

workflow_results <- list(
    mode = workflow_mode,
    easy = easy_results,
    hard = hard_results,
    hard_sampling_approved = hard_sampling_approved
)
