source("PDBEntryBuilder_v2.R")

entry <- PDBEntryBuilder$new(
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
            dimensions = list(beta = 2, sigma = 1)
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
            dimensions = list(beta = 2, sigma = 1)
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
    continue_on_error = TRUE
)

entry$summarize_workflow_results(results)
