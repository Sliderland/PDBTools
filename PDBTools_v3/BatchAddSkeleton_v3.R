.pdbtools_v3_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.pdbtools_v3_file) || !nzchar(.pdbtools_v3_file)) {
  .pdbtools_v3_file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(.pdbtools_v3_file_arg)) {
    .pdbtools_v3_file <- sub("^--file=", "", .pdbtools_v3_file_arg[[1L]])
  }
}
if (is.null(.pdbtools_v3_file) || !nzchar(.pdbtools_v3_file)) {
  stop("Run or source this workflow from its file path.", call. = FALSE)
}
.pdbtools_v3_dir <- dirname(normalizePath(.pdbtools_v3_file, mustWork = TRUE))
.pdbtools_repo_root <- normalizePath(
  file.path(.pdbtools_v3_dir, ".."),
  mustWork = TRUE
)
source(file.path(.pdbtools_v3_dir, "PDBEntryBuilder_v3.R"))

entry <- PDBEntryBuilder_v3$new(
    "/Users/gerpr308/Documents/posteriordb"
)

common_model <- list(
    stan_file = source_model_file,
    info = test_model_info
)

entry_specs <- list(
    earnings_1 = list(
        data = list(
            data = earnings_data,
            info = make_data_info("pdbtools_earnings_1", 1)
        ),
        posterior = list(
            parameters_include = c("beta", "sigma")
        ),
        reference = list(
            sampling_args = sampling_args,
            comments = "First integration-test posterior."
        )
    ),
    earnings_2 = list(
        data = list(
            data = earnings_data,
            info = make_data_info("pdbtools_earnings_2", 2)
        ),
        posterior = list(
            parameters_include = c("beta", "sigma")
        ),
        reference = list(
            sampling_args = sampling_args,
            comments = "Second integration-test posterior."
        )
    )
)

results <- entry$run_workflows(
    model = common_model,
    entries = entry_specs,
    register = TRUE,
    sample = TRUE,
    write = TRUE,
    overwrite = TRUE,
    continue_on_error = TRUE,
    save_failed_fits = FALSE,
    save_all_fits = FALSE,
    failed_fit_dir = path.expand(
        "~/Documents/PDBTools_diagnostics/failed_reference_fits"
    ),
    all_fit_dir = path.expand(
        "~/Documents/PDBTools_diagnostics/all_reference_fits"
    )
)

entry$summarize_workflow_results(results)
