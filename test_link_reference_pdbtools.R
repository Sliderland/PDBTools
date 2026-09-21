# Regression test for linking reference files to an existing posterior.
# Run from this directory with:
#   Rscript test_link_reference_pdbtools.R

source("PDBEntryBuilder_v2.R")

local({
    root <- tempfile("pdbtools-link-")
    database <- file.path(root, "posterior_database")
    dirs <- c(
        "data/info",
        "data/data",
        "models/info",
        "models/stan",
        "posteriors",
        "reference_posteriors/draws/info",
        "reference_posteriors/draws/draws",
        "alias"
    )
    dir.create(database, recursive = TRUE)
    for (directory in dirs) {
        dir.create(file.path(database, directory), recursive = TRUE)
    }
    writeLines("{}", file.path(database, "alias/posteriors.json"))

    write_json <- function(object, path) {
        jsonlite::write_json(
            object,
            path,
            pretty = TRUE,
            auto_unbox = TRUE,
            null = "null",
            digits = NA
        )
    }
    write_zip <- function(directory, filename, contents) {
        archive <- file.path(database, directory, filename)
        source <- sub("[.]zip$", "", archive)
        writeLines(contents, source)
        oldwd <- setwd(dirname(archive))
        on.exit(setwd(oldwd), add = TRUE)
        utils::zip(basename(archive), basename(source), flags = "-jq")
        setwd(oldwd)
        unlink(source)
    }

    write_json(
        list(
            name = "toy_data",
            title = "Toy data",
            data_file = "data/data/toy_data.json",
            added_by = "test",
            added_date = "2026-01-01"
        ),
        file.path(database, "data/info/toy_data.info.json")
    )
    write_zip("data/data", "toy_data.json.zip", '{"y":[1]}')

    write_json(
        list(
            name = "toy_model",
            title = "Toy model",
            added_by = "test",
            added_date = "2026-01-01",
            model_implementations = list(
                stan = list(model_code = "models/stan/toy_model.stan")
            )
        ),
        file.path(database, "models/info/toy_model.info.json")
    )
    writeLines("parameters {}", file.path(database, "models/stan/toy_model.stan"))

    posterior_name <- "toy_data-toy_model"
    write_json(
        list(
            name = posterior_name,
            model_name = "toy_model",
            data_name = "toy_data",
            reference_posterior_name = NULL,
            dimensions = list(y = 1),
            added_by = "test",
            added_date = "2026-01-01"
        ),
        file.path(database, "posteriors", paste0(posterior_name, ".json"))
    )
    write_json(
        list(
            name = posterior_name,
            inference = list(method = "stan_sampling", method_arguments = list()),
            diagnostics = list(ndraws = 1),
            checks_made = list(),
            comments = NULL,
            added_by = "test",
            added_date = "2026-01-01",
            versions = list()
        ),
        file.path(
            database,
            "reference_posteriors/draws/info",
            paste0(posterior_name, ".info.json")
        )
    )
    write_zip(
        "reference_posteriors/draws/draws",
        paste0(posterior_name, ".json.zip"),
        '[{"y":1}]'
    )

    builder <- PDBEntryBuilder$new(root, detect_cores = FALSE, n_cores = 1)
    changed <- builder$link_reference_posterior(posterior_name)
    stopifnot(isTRUE(changed))

    posterior_after <- jsonlite::read_json(
        file.path(database, "posteriors", paste0(posterior_name, ".json")),
        simplifyVector = FALSE
    )
    stopifnot(
        identical(
            posterior_after$reference_posterior_name,
            posterior_name
        )
    )

    changed_again <- builder$link_reference_posterior(posterior_name)
    stopifnot(identical(changed_again, FALSE))

    unlink(root, recursive = TRUE, force = TRUE)
    message("PDBTools reference-link regression test passed.")
})
