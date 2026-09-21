# Focused regression test for PDBEntryBuilder$rename_pdb().
# Run from this directory with:
#   Rscript test_rename_pdbtools.R

source("PDBEntryBuilder_v2.R")

local({
    root <- tempfile("pdbtools-rename-")
    database <- file.path(root, "posterior_database")
    dirs <- c(
        "data/info",
        "data/data",
        "models/info",
        "models/stan",
        "models/r",
        "posteriors",
        "reference_posteriors/draws/info",
        "reference_posteriors/draws/draws",
        "reference_posteriors/summary_statistics/mean_value/info",
        "reference_posteriors/summary_statistics/mean_value/mean_value",
        "reference_posteriors/summary_statistics/mean_squared_value/info",
        "reference_posteriors/summary_statistics/mean_squared_value/mean_squared_value",
        "alias"
    )
    dir.create(database, recursive = TRUE)
    for (directory in dirs) {
        dir.create(file.path(database, directory), recursive = TRUE)
    }

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

    data_info <- list(
        name = "data_old",
        title = "Data",
        added_by = "test",
        added_date = "2026-01-01",
        data_file = "data/data/data_old.json"
    )
    write_json(data_info, file.path(database, "data/info/data_old.info.json"))
    write_zip("data/data", "data_old.json.zip", '{"y":[1,2]}')
    writeLines('{"y":[1,2]}', file.path(database, "data/data/data_old.json"))

    model_info <- list(
        name = "model_old",
        title = "Model",
        added_by = "test",
        added_date = "2026-01-01",
        model_implementations = list(
            stan = list(model_code = "models/stan/model_old.stan"),
            r = list(model_code = "models/r/model_old.R")
        )
    )
    write_json(model_info, file.path(database, "models/info/model_old.info.json"))
    writeLines("parameters {}", file.path(database, "models/stan/model_old.stan"))
    writeLines("model_old <- function(x) x", file.path(database, "models/r/model_old.R"))

    posterior_info <- list(
        name = "data_old-model_old",
        model_name = "model_old",
        data_name = "data_old",
        reference_posterior_name = "data_old-model_old",
        dimensions = list(y = 2),
        added_by = "test",
        added_date = "2026-01-01"
    )
    write_json(
        posterior_info,
        file.path(database, "posteriors/data_old-model_old.json")
    )

    reference_info <- list(
        name = "data_old-model_old",
        inference = list(method = "analytical", method_arguments = list()),
        diagnostics = list(ndraws = 1),
        checks_made = list(),
        comments = "test",
        added_by = "test",
        added_date = "2026-01-01",
        versions = list()
    )
    write_json(
        reference_info,
        file.path(
            database,
            "reference_posteriors/draws/info/data_old-model_old.info.json"
        )
    )
    write_zip(
        "reference_posteriors/draws/draws",
        "data_old-model_old.json.zip",
        '[{"y":[1]}]'
    )
    for (summary_type in c("mean_value", "mean_squared_value")) {
        write_json(
            reference_info,
            file.path(
                database,
                "reference_posteriors",
                "summary_statistics",
                summary_type,
                "info",
                "data_old-model_old.info.json"
            )
        )
        write_json(
            list(names = "y", value = 1),
            file.path(
                database,
                "reference_posteriors",
                "summary_statistics",
                summary_type,
                summary_type,
                "data_old-model_old.json"
            )
        )
    }
    write_json(
        list(
            "data_old-model_old" = "data_old-model_old",
            alias_for_posterior = "data_old-model_old"
        ),
        file.path(database, "alias/posteriors.json")
    )

    builder <- PDBEntryBuilder$new(root, detect_cores = FALSE, n_cores = 1)

    builder$rename_pdb("data_old", "data_new", type = "data")
    stopifnot(
        file.exists(file.path(database, "data/info/data_new.info.json")),
        file.exists(file.path(database, "data/data/data_new.json")),
        file.exists(file.path(database, "data/data/data_new.json.zip")),
        !file.exists(file.path(database, "data/info/data_old.info.json")),
        !file.exists(file.path(database, "data/data/data_old.json")),
        !file.exists(file.path(database, "data/data/data_old.json.zip"))
    )
    stopifnot(
        jsonlite::read_json(
            file.path(database, "data/info/data_new.info.json"),
            simplifyVector = TRUE
        )$data_file == "data/data/data_new.json",
        utils::unzip(
            file.path(database, "data/data/data_new.json.zip"),
            list = TRUE
        )$Name == "data_new.json"
    )

    builder$rename_pdb("model_old", "model_new", type = "model")
    stopifnot(
        file.exists(file.path(database, "models/stan/model_new.stan")),
        file.exists(file.path(database, "models/r/model_new.R")),
        !file.exists(file.path(database, "models/stan/model_old.stan")),
        !file.exists(file.path(database, "models/r/model_old.R"))
    )
    model_after <- jsonlite::read_json(
        file.path(database, "models/info/model_new.info.json"),
        simplifyVector = FALSE
    )
    stopifnot(
        model_after$name == "model_new",
        model_after$model_implementations$stan$model_code ==
            "models/stan/model_new.stan",
        model_after$model_implementations$r$model_code ==
            "models/r/model_new.R"
    )

    posterior_after <- jsonlite::read_json(
        file.path(database, "posteriors/data_new-model_new.json"),
        simplifyVector = FALSE
    )
    stopifnot(
        posterior_after$name == "data_new-model_new",
        posterior_after$data_name == "data_new",
        posterior_after$model_name == "model_new",
        posterior_after$reference_posterior_name == "data_new-model_new",
        file.exists(file.path(
            database,
            "reference_posteriors/draws/info/data_new-model_new.info.json"
        )),
        file.exists(file.path(
            database,
            "reference_posteriors/summary_statistics/mean_value/mean_value/data_new-model_new.json"
        )),
        file.exists(file.path(
            database,
            "reference_posteriors/summary_statistics/mean_squared_value/mean_squared_value/data_new-model_new.json"
        )),
        utils::unzip(
            file.path(
                database,
                "reference_posteriors/draws/draws/data_new-model_new.json.zip"
            ),
            list = TRUE
        )$Name == "data_new-model_new.json"
    )
    aliases <- jsonlite::read_json(
        file.path(database, "alias/posteriors.json"),
        simplifyVector = TRUE
    )
    stopifnot(
        aliases[["data_old-model_old"]] == "data_new-model_new",
        aliases[["alias_for_posterior"]] == "data_new-model_new"
    )

    unchanged <- readLines(file.path(database, "posteriors/data_new-model_new.json"))
    failed <- try(builder$rename_pdb("data_new", "../bad", type = "data"), silent = TRUE)
    stopifnot(inherits(failed, "try-error"))
    stopifnot(file.exists(file.path(database, "posteriors/data_new-model_new.json")))
    stopifnot(identical(readLines(file.path(database, "posteriors/data_new-model_new.json")), unchanged))

    write_json(
        list(name = "data_taken"),
        file.path(database, "data/info/data_taken.info.json")
    )
    failed_collision <- try(
        builder$rename_pdb("data_new", "data_taken", type = "data"),
        silent = TRUE
    )
    stopifnot(inherits(failed_collision, "try-error"))
    stopifnot(file.exists(file.path(database, "data/info/data_new.info.json")))

    unlink(root, recursive = TRUE, force = TRUE)
    message("PDBTools rename regression test passed.")
})
