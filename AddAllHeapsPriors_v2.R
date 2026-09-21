# Register and sample the models distributed in HeapsStanPrograms.
#
# Models:
# - statprior: stationary VAR with the exchangeable prior from Heaps (2023)
# - semiconj: non-stationary semi-conjugate comparison model
# - statrmlprior: stationary VAR using the Roy et al. reparameterisation
# - statinvertprior: stationary and invertible VARMA(p, q)
# - ansleykohn: stationary VAR using the Ansley-Kohn reparameterisation
#
# Each model family has a distinct Stan data block. Consequently, this script
# creates model-specific PosteriorDB data entries rather than passing one
# oversized or incorrectly shaped data list to every model.

source("PDBEntryBuilder_v2.R")

pdb_path <- path.expand("~/Documents/posteriordb")
heaps_program_path <- normalizePath(
  file.path("HeapsStanPrograms"),
  mustWork = TRUE
)

# Execution controls -------------------------------------------------------

register_entries <- TRUE
run_sampling <- TRUE
write_reference_files <- TRUE

overwrite_registration <- TRUE
overwrite_reference_files <- FALSE
skip_completed_references <- TRUE
continue_on_error <- TRUE

# Retain completed but non-passing rstanfit objects outside PosteriorDB for
# later diagnostics. These files can be large; leave disabled unless needed.
save_failed_fits <- FALSE
failed_fit_dir <- path.expand(
  "~/Documents/PDBTools_diagnostics/failed_reference_fits"
)

# Run the small and medium-dimensional datasets. Add "med20" back for the
# substantially longer large-dimensional jobs.
dataset_sizes_to_run <- c("small3", "med10")
models_to_run <- c(
  "statprior_var",
  "semiconj_var",
  "statrml_var",
  "statinvert_varma",
  "ansleykohn_var"
)

if (write_reference_files && !run_sampling) {
  stop("`write_reference_files = TRUE` requires `run_sampling = TRUE`.")
}

entry <- PDBEntryBuilder$new(pdb_path)

# Bibliography -------------------------------------------------------------

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

# Observed time-series data ------------------------------------------------

source(file.path(heaps_program_path, "read.R"))

dataset_definitions <- list(
  small3 = list(dimension = 3L, size = "Small"),
  med10 = list(dimension = 10L, size = "Medium"),
  med20 = list(dimension = 20L, size = "Large")
)

unknown_datasets <- setdiff(
  dataset_sizes_to_run,
  names(dataset_definitions)
)
if (length(unknown_datasets) > 0L) {
  stop(
    "Unknown dataset selections: ",
    paste(unknown_datasets, collapse = ", "),
    call. = FALSE
  )
}
dataset_definitions <- dataset_definitions[dataset_sizes_to_run]

build_observations <- function(dimension) {
  process_data(
    num = dimension,
    omit = c(1:2, 197:200),
    Nahead = 40,
    data_dir_path = file.path(heaps_program_path, "data")
  )$y
}

observations <- lapply(dataset_definitions, function(definition) {
  build_observations(definition$dimension)
})

# Prior hyperparameters ----------------------------------------------------

# These values reproduce heaps_list() in helper.R, bridgestanHelper.R, and
# generateComparisonList.R in the Bayesian Time Series project.
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

base_stan_data <- function(y) {
  list(
    m = ncol(y),
    p = prior_defaults$p,
    N = nrow(y),
    y = y
  )
}

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

validate_model_data <- function(model_name, data) {
  expected <- expected_data_fields[[model_name]]
  if (is.null(expected)) {
    stop("Unknown model: ", model_name, call. = FALSE)
  }

  missing_fields <- setdiff(expected, names(data))
  unnecessary_fields <- setdiff(names(data), expected)
  if (length(missing_fields) > 0L || length(unnecessary_fields) > 0L) {
    problems <- c(
      if (length(missing_fields) > 0L) {
        paste0("missing: ", paste(missing_fields, collapse = ", "))
      },
      if (length(unnecessary_fields) > 0L) {
        paste0(
          "unnecessary: ",
          paste(unnecessary_fields, collapse = ", ")
        )
      }
    )
    stop(
      "Invalid data list for `", model_name, "` (",
      paste(problems, collapse = "; "), ").",
      call. = FALSE
    )
  }

  data
}

build_model_data <- function(model_name, y) {
  base <- base_stan_data(y)
  covariance_prior <- list(
    scale_diag = prior_defaults$scale_diag,
    scale_offdiag = prior_defaults$scale_offdiag,
    df = ncol(y) + 4
  )

  data <- switch(model_name,
    statprior_var = c(
      base,
      list(
        es = prior_defaults$es,
        fs = prior_defaults$fs,
        gs = prior_defaults$gs,
        hs = prior_defaults$hs
      ),
      covariance_prior
    ),
    semiconj_var = c(
      base,
      list(
        m_diag = prior_defaults$m_diag,
        s_diag = prior_defaults$s_diag,
        m_offdiag = prior_defaults$m_offdiag,
        s_offdiag = prior_defaults$s_offdiag
      ),
      covariance_prior
    ),
    statrml_var = c(
      base,
      list(df = ncol(y) + 4)
    ),
    statinvert_varma = c(
      base,
      list(
        q = prior_defaults$q,
        # array[2] vector[2] in Stan: rows correspond to the AR and
        # MA components; columns correspond to diagonal/off-diagonal.
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

  validate_model_data(model_name, data)
}

# Model definitions --------------------------------------------------------

generated_quantity_exclusions <- c(
  "topblock",
  "companion",
  "lambdas",
  "lambda_moduli",
  "max_lambda_modulus"
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
    keywords = c(
      "stationary", "invertible", "VARMA", "exchangeable"
    )
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

unknown_models <- setdiff(models_to_run, names(model_definitions))
if (length(unknown_models) > 0L) {
  stop(
    "Unknown model selections: ",
    paste(unknown_models, collapse = ", "),
    call. = FALSE
  )
}
model_definitions <- model_definitions[models_to_run]

model_info <- Map(
  function(name, definition) {
    list(
      name = name,
      keywords = definition$keywords,
      title = definition$title,
      description = definition$description,
      urls = "https://doi.org/10.1080/10618600.2022.2079648",
      references = "heaps2023stationary",
      framework = "stan",
      added_by = "Gerald Press",
      added_date = Sys.Date()
    )
  },
  names(model_definitions),
  model_definitions
)

# Create model-specific data entries --------------------------------------

data_entries <- list()
for (dataset_key in names(dataset_definitions)) {
  definition <- dataset_definitions[[dataset_key]]
  y <- observations[[dataset_key]]

  for (model_name in names(model_definitions)) {
    # Preserve the established data names for statprior_var. Other model
    # families receive suffixes because their hyperparameter fields differ.
    data_name <- if (identical(model_name, "statprior_var")) {
      paste0("heaps_", dataset_key)
    } else {
      paste0("heaps_", dataset_key, "_", model_name)
    }

    data_entries[[data_name]] <- list(
      model_name = model_name,
      data = build_model_data(model_name, y),
      info = list(
        name = data_name,
        keywords = c(
          "United States",
          "Macroeconomics",
          "Vector Autoregression",
          definition$size,
          model_name
        ),
        title = paste(
          definition$size,
          "Heaps macroeconomic data for",
          model_name
        ),
        description = paste0(
          "Quarterly US macroeconomic observations with ",
          definition$dimension,
          " variables and the hyperparameters required by ",
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
    )
  }
}

plan <- data.frame(
  data_name = names(data_entries),
  model_name = vapply(
    data_entries,
    function(x) x$model_name,
    character(1)
  ),
  stringsAsFactors = FALSE
)
plan$posterior_name <- paste(plan$data_name, plan$model_name, sep = "-")
print(plan)

# Registration -------------------------------------------------------------

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

exists_exactly <- function(matches, name) name %in% matches

if (register_entries) {
  registration_results$bibliography <- attempt(
    "Bibliography registration",
    entry$add_bibtex_entry(heaps_reference)
  )

  for (data_name in names(data_entries)) {
    existing <- exists_exactly(entry$search_data(data_name), data_name)
    registration_results$data[[data_name]] <- if (
      existing && !overwrite_registration
    ) {
      list(status = "skipped_existing", value = NULL, error = NULL)
    } else {
      attempt(
        paste("Data registration", data_name),
        entry$add_data(
          data = data_entries[[data_name]]$data,
          info = data_entries[[data_name]]$info,
          overwrite = overwrite_registration
        )
      )
    }
  }

  for (model_name in names(model_definitions)) {
    existing <- exists_exactly(entry$search_model(model_name), model_name)
    registration_results$models[[model_name]] <- if (
      existing && !overwrite_registration
    ) {
      list(status = "skipped_existing", value = NULL, error = NULL)
    } else {
      attempt(
        paste("Model registration", model_name),
        entry$add_model_code(
          stan_file = model_definitions[[model_name]]$stan_file,
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
    existing <- exists_exactly(
      entry$search_posterior(posterior_name),
      posterior_name
    )

    registration_results$posteriors[[posterior_name]] <- if (
      existing && !overwrite_registration
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
            parameters_exclude = generated_quantity_exclusions
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

# Sequential reference sampling ------------------------------------------

sampling_args <- list(
  chains = 10L,
  iter = 30000L,
  warmup = 10000L,
  refresh = 10000L,
  thin = 20L,
  seed = 123L,
  control = list(adapt_delta = 0.95)
)

sampling_specs <- setNames(
  lapply(seq_len(nrow(plan)), function(row) {
    list(
      posterior_name = plan$posterior_name[[row]],
      posterior = list(
        data_name = plan$data_name[[row]],
        model_name = plan$model_name[[row]]
      ),
      reference = list(
        sampling_args = sampling_args,
        comments = paste(
          "Sequential reference sampling for",
          plan$posterior_name[[row]]
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
    "Sampling is disabled. Review `plan`, then enable registration and ",
    "sampling in separate runs."
  )
}
