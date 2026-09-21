# The PosteriorDB name is part of the database graph, not just a filename.
# These helpers stage a complete rename before changing any source files.

pdbtools_rename_validate_name <- function(name) {
    if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
        stop("Names must be non-empty character strings.", call. = FALSE)
    }
    if (
        grepl("[/\\\\]", name) ||
            name %in% c(".", "..") ||
            grepl("[[:cntrl:]]", name)
    ) {
        stop(
            "Names must be single path components without separators or control characters.",
            call. = FALSE
        )
    }
    invisible(name)
}

pdbtools_rename_relpath <- function(root, path) {
    path <- gsub("\\\\", "/", path)
    if (
        !nzchar(path) ||
            grepl("^[/]", path) ||
            grepl("(^|/)[.][.](/|$)", path)
    ) {
        stop("PosteriorDB paths must be relative and cannot contain '..'.", call. = FALSE)
    }
    file.path(root, path)
}

pdbtools_rename_read_json <- function(path) {
    tryCatch(
        jsonlite::read_json(path, simplifyVector = FALSE),
        error = function(e) {
            stop(
                "Cannot parse JSON file '",
                path,
                "': ",
                conditionMessage(e),
                call. = FALSE
            )
        }
    )
}

pdbtools_rename_json_text <- function(object) {
    out <- jsonlite::toJSON(
        object,
        pretty = TRUE,
        auto_unbox = TRUE,
        null = "null",
        digits = NA,
        encoding = "UTF-8"
    )
    Encoding(out) <- "UTF-8"
    as.character(out)
}

pdbtools_rename_model_path <- function(value, old_name, new_name) {
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
        return(value)
    }
    value <- gsub("\\\\", "/", value)
    if (!startsWith(value, "models/") ||
        !startsWith(basename(value), paste0(old_name, "."))) {
        return(value)
    }
    paste0(
        dirname(value),
        "/",
        new_name,
        substring(basename(value), nchar(old_name) + 1L)
    )
}

pdbtools_rename_add_action <- function(
    actions,
    source,
    target,
    content = NULL,
    zip_member = NULL
) {
    source <- gsub("\\\\", "/", source)
    target <- gsub("\\\\", "/", target)
    previous <- actions[[source]]
    if (!is.null(previous)) {
        same_content <-
            (is.null(previous$content) && is.null(content)) ||
                identical(previous$content, content)
        if (
            !identical(previous$target, target) ||
                !same_content ||
                !identical(previous$zip_member, zip_member)
        ) {
            stop(
                "The rename plan contains conflicting actions for '",
                source,
                "'.",
                call. = FALSE
            )
        }
        return(actions)
    }
    actions[[source]] <- list(
        source = source,
        target = target,
        content = content,
        zip_member = zip_member
    )
    actions
}

pdbtools_rename_stage_zip <- function(source, destination, member_name) {
    listing <- tryCatch(
        utils::unzip(source, list = TRUE)$Name,
        error = function(e) {
            stop(
                "Cannot inspect ZIP archive '",
                source,
                "': ",
                conditionMessage(e),
                call. = FALSE
            )
        }
    )
    files <- listing[!grepl("/$", listing)]
    if (length(files) != 1L || !grepl("[.]json$", files[[1L]])) {
        stop(
            "Expected a ZIP archive containing exactly one JSON file: ",
            source,
            call. = FALSE
        )
    }

    extraction <- tempfile(".pdbtools-rename-unzip-")
    dir.create(extraction, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(extraction, recursive = TRUE, force = TRUE), add = TRUE)
    utils::unzip(source, exdir = extraction)
    old_member <- file.path(extraction, files[[1L]])
    new_member <- file.path(extraction, member_name)
    dir.create(dirname(new_member), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(old_member, new_member)) {
        stop("Could not rename the ZIP member.", call. = FALSE)
    }

    oldwd <- setwd(extraction)
    on.exit(setwd(oldwd), add = TRUE)
    zip_destination <- paste0(destination, ".zip")
    utils::zip(zip_destination, files = member_name, flags = "-jq")
    setwd(oldwd)
    if (!file.rename(zip_destination, destination)) {
        stop("Could not finalize the staged ZIP archive.", call. = FALSE)
    }
    invisible(TRUE)
}

pdbtools_rename_commit <- function(actions, root) {
    if (!length(actions)) {
        return(invisible(TRUE))
    }

    action_list <- unname(actions)
    sources <- vapply(action_list, `[[`, character(1), "source")
    targets <- vapply(action_list, `[[`, character(1), "target")
    if (anyDuplicated(targets)) {
        stop("The rename plan has duplicate target paths.", call. = FALSE)
    }

    source_abs <- vapply(
        sources,
        function(path) pdbtools_rename_relpath(root, path),
        character(1)
    )
    target_abs <- vapply(
        targets,
        function(path) pdbtools_rename_relpath(root, path),
        character(1)
    )
    if (any(!file.exists(source_abs))) {
        missing <- sources[!file.exists(source_abs)]
        stop(
            "The rename plan refers to missing source file(s): ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    source_set <- unique(source_abs)
    unexpected_targets <- target_abs[
        file.exists(target_abs) & !(target_abs %in% source_set)
    ]
    if (length(unexpected_targets)) {
        stop(
            "Refusing to overwrite existing file '",
            unexpected_targets[[1L]],
            "'.",
            call. = FALSE
        )
    }
    if (any(target_abs %in% source_abs & target_abs != source_abs)) {
        stop(
            "The rename plan contains a target that is another source file.",
            call. = FALSE
        )
    }

    stage <- tempfile(".pdbtools-rename-stage-", tmpdir = dirname(root))
    backup <- tempfile(".pdbtools-rename-backup-", tmpdir = dirname(root))
    dir.create(stage, recursive = TRUE, showWarnings = FALSE)
    dir.create(backup, recursive = TRUE, showWarnings = FALSE)
    cleanup <- function() {
        unlink(stage, recursive = TRUE, force = TRUE)
        unlink(backup, recursive = TRUE, force = TRUE)
    }
    on.exit(cleanup(), add = TRUE)

    staged <- character(length(action_list))
    backups <- character(length(action_list))
    for (i in seq_along(action_list)) {
        staged[[i]] <- file.path(stage, sprintf("%06d", i))
        backups[[i]] <- file.path(backup, sprintf("%06d", i))
        action <- action_list[[i]]
        if (!is.null(action$zip_member)) {
            pdbtools_rename_stage_zip(
                source_abs[[i]],
                staged[[i]],
                action$zip_member
            )
        } else if (!is.null(action$content)) {
            writeLines(action$content, staged[[i]], useBytes = TRUE)
        } else if (!file.copy(source_abs[[i]], staged[[i]], overwrite = FALSE)) {
            stop("Could not stage '", sources[[i]], "'.", call. = FALSE)
        }
    }

    backed_up <- 0L
    committed <- 0L
    rollback <- function() {
        if (committed > 0L) {
            for (i in seq.int(committed, 1L)) {
                unlink(target_abs[[i]], force = TRUE)
            }
        }
        if (backed_up > 0L) {
            for (i in seq.int(backed_up, 1L)) {
                if (file.exists(backups[[i]])) {
                    dir.create(
                        dirname(source_abs[[i]]),
                        recursive = TRUE,
                        showWarnings = FALSE
                    )
                    file.rename(backups[[i]], source_abs[[i]])
                }
            }
        }
    }

    for (i in seq_along(action_list)) {
        if (!file.rename(source_abs[[i]], backups[[i]])) {
            rollback()
            stop(
                "Could not reserve source file '",
                sources[[i]],
                "'; no changes were kept.",
                call. = FALSE
            )
        }
        backed_up <- i
    }
    for (i in seq_along(action_list)) {
        dir.create(
            dirname(target_abs[[i]]),
            recursive = TRUE,
            showWarnings = FALSE
        )
        if (!file.rename(staged[[i]], target_abs[[i]])) {
            rollback()
            stop(
                "Could not install target file '",
                targets[[i]],
                "'; the migration was rolled back.",
                call. = FALSE
            )
        }
        committed <- i
    }
    invisible(TRUE)
}

pdbtools_rename_reference_files <- function(root, old_name) {
    reference_root <- file.path(root, "reference_posteriors")
    if (!dir.exists(reference_root)) {
        return(character())
    }
    files <- list.files(reference_root, full.names = TRUE, recursive = TRUE)
    files[basename(files) %in% c(
        paste0(old_name, ".info.json"),
        paste0(old_name, ".json"),
        paste0(old_name, ".json.zip")
    )]
}

pdbtools_rename_add_reference_actions <- function(
    actions,
    root,
    old_name,
    new_name
) {
    files <- pdbtools_rename_reference_files(root, old_name)
    if (!length(files)) {
        stop(
            "Reference posterior '",
            old_name,
            "' is linked but no reference files were found.",
            call. = FALSE
        )
    }
    for (source in files) {
        relative <- substring(
            normalizePath(source, mustWork = FALSE),
            nchar(normalizePath(root, mustWork = FALSE)) + 2L
        )
        target <- file.path(
            dirname(relative),
            paste0(new_name, substring(basename(relative), nchar(old_name) + 1L))
        )
        content <- NULL
        if (grepl("[.]info[.]json$", relative)) {
            info <- pdbtools_rename_read_json(source)
            info$name <- new_name
            content <- pdbtools_rename_json_text(info)
        }
        zip_member <- if (grepl("[.]json[.]zip$", relative)) {
            paste0(new_name, ".json")
        } else {
            NULL
        }
        actions <- pdbtools_rename_add_action(
            actions,
            relative,
            target,
            content,
            zip_member
        )
    }
    actions
}

pdbtools_rename_infer_type <- function(root, old_name, type = NULL) {
    choices <- c("data", "model", "posterior")
    if (!is.null(type)) {
        if (!is.character(type) || length(type) != 1L ||
            is.na(type) || !type %in% choices) {
            stop("`type` must be one of: data, model, posterior.", call. = FALSE)
        }
        return(type)
    }
    exists <- c(
        data = file.exists(file.path(root, "data", "info", paste0(old_name, ".info.json"))),
        model = file.exists(file.path(root, "models", "info", paste0(old_name, ".info.json"))),
        posterior = file.exists(file.path(root, "posteriors", paste0(old_name, ".json")))
    )
    if (sum(exists) != 1L) {
        stop(
            "`type` is required unless the old name identifies exactly one data, model, or posterior.",
            call. = FALSE
        )
    }
    names(exists)[exists]
}

pdbtools_rename_entity <- function(
    database_path,
    old_name,
    new_name,
    type = NULL
) {
    root <- normalizePath(
        file.path(database_path, "posterior_database"),
        mustWork = TRUE
    )
    pdbtools_rename_validate_name(old_name)
    pdbtools_rename_validate_name(new_name)
    if (identical(old_name, new_name)) {
        return(list(
            type = pdbtools_rename_infer_type(root, old_name, type),
            old_name = old_name,
            new_name = new_name,
            posterior_renames = list(),
            affected_posteriors = character(),
            files = character()
        ))
    }
    type <- pdbtools_rename_infer_type(root, old_name, type)
    actions <- list()

    add_action <- function(source, target, content = NULL, zip_member = NULL) {
        actions <<- pdbtools_rename_add_action(
            actions,
            source,
            target,
            content,
            zip_member
        )
    }

    if (identical(type, "data")) {
        info_rel <- file.path("data", "info", paste0(old_name, ".info.json"))
        info_path <- pdbtools_rename_relpath(root, info_rel)
        data_rel <- file.path(
            "data",
            "data",
            paste0(old_name, c(".json.zip", ".json"))
        )
        data_rel <- data_rel[vapply(
            data_rel,
            function(path) file.exists(pdbtools_rename_relpath(root, path)),
            logical(1)
        )]
        if (!length(data_rel)) {
            stop("No data file was found for '", old_name, "'.", call. = FALSE)
        }
        data_info <- pdbtools_rename_read_json(info_path)
        if (!identical(data_info$name, old_name)) {
            stop("Data metadata name does not match its filename.", call. = FALSE)
        }
        data_info$name <- new_name
        # PosteriorDB metadata uses the unzipped logical JSON path even when
        # the stored data object is compressed as <name>.json.zip.
        data_info$data_file <- paste0("data/data/", new_name, ".json")
        add_action(
            info_rel,
            file.path("data", "info", paste0(new_name, ".info.json")),
            pdbtools_rename_json_text(data_info)
        )
        for (source in data_rel) {
            target <- file.path(
                dirname(source),
                paste0(new_name, substring(basename(source), nchar(old_name) + 1L))
            )
            member <- if (grepl("[.]json[.]zip$", source)) {
                paste0(new_name, ".json")
            } else {
                NULL
            }
            add_action(source, target, zip_member = member)
        }
        entity_field <- "data_name"
    } else if (identical(type, "model")) {
        info_rel <- file.path("models", "info", paste0(old_name, ".info.json"))
        info_path <- pdbtools_rename_relpath(root, info_rel)
        model_info <- pdbtools_rename_read_json(info_path)
        if (!identical(model_info$name, old_name)) {
            stop("Model metadata name does not match its filename.", call. = FALSE)
        }
        implementations <- model_info$model_implementations
        if (is.null(implementations) || !length(implementations)) {
            stop("Model '", old_name, "' has no implementations.", call. = FALSE)
        }
        model_info$name <- new_name
        for (implementation in seq_along(implementations)) {
            fields <- names(implementations[[implementation]])
            for (field in fields) {
                value <- implementations[[implementation]][[field]]
                new_value <- pdbtools_rename_model_path(value, old_name, new_name)
                if (!identical(value, new_value)) {
                    source <- gsub("\\\\", "/", value)
                    target <- gsub("\\\\", "/", new_value)
                    pdbtools_rename_relpath(root, source)
                    pdbtools_rename_relpath(root, target)
                    add_action(source, target)
                    model_info$model_implementations[[implementation]][[field]] <- new_value
                }
            }
        }
        add_action(
            info_rel,
            file.path("models", "info", paste0(new_name, ".info.json")),
            pdbtools_rename_json_text(model_info)
        )
        entity_field <- "model_name"
    } else {
        posterior_rel <- file.path("posteriors", paste0(old_name, ".json"))
        if (!file.exists(pdbtools_rename_relpath(root, posterior_rel))) {
            stop("No posterior was found for '", old_name, "'.", call. = FALSE)
        }
        entity_field <- NULL
    }

    posterior_dir <- file.path(root, "posteriors")
    posterior_paths <- if (dir.exists(posterior_dir)) {
        list.files(
            posterior_dir,
            pattern = "[.]json$",
            full.names = TRUE,
            recursive = FALSE
        )
    } else {
        character()
    }
    posterior_records <- lapply(posterior_paths, function(path) {
        record <- pdbtools_rename_read_json(path)
        if (is.null(record$name)) {
            stop("Posterior file '", path, "' has no `name` field.", call. = FALSE)
        }
        record$.rename_path <- path
        record
    })

    affected <- list()
    for (record in posterior_records) {
        is_affected <- if (is.null(entity_field)) {
            identical(record$name, old_name)
        } else {
            identical(record[[entity_field]], old_name)
        }
        if (!is_affected) {
            next
        }

        old_posterior_name <- record$name
        other_name <- if (identical(type, "data")) {
            record$model_name
        } else {
            record$data_name
        }
        expected_old <- if (is.null(entity_field)) {
            old_name
        } else if (identical(type, "data")) {
            paste(old_name, other_name, sep = "-")
        } else {
            paste(other_name, old_name, sep = "-")
        }
        new_posterior_name <- old_posterior_name
        if (identical(old_posterior_name, expected_old)) {
            new_posterior_name <- if (identical(type, "data")) {
                paste(new_name, other_name, sep = "-")
            } else if (identical(type, "model")) {
                paste(other_name, new_name, sep = "-")
            } else {
                new_name
            }
        }
        if (identical(type, "data")) {
            record$data_name <- new_name
        }
        if (identical(type, "model")) {
            record$model_name <- new_name
        }
        if (
            identical(record$reference_posterior_name, old_posterior_name) &&
                !identical(new_posterior_name, old_posterior_name)
        ) {
            record$reference_posterior_name <- new_posterior_name
        }
        record$name <- new_posterior_name

        actual_relative <- substring(
            normalizePath(record$.rename_path, mustWork = FALSE),
            nchar(normalizePath(root, mustWork = FALSE)) + 2L
        )
        expected_relative <- file.path(
            "posteriors",
            paste0(old_posterior_name, ".json")
        )
        if (!identical(actual_relative, expected_relative)) {
            stop(
                "Posterior file '",
                actual_relative,
                "' disagrees with its `name` field.",
                call. = FALSE
            )
        }
        serializable <- record[setdiff(names(record), ".rename_path")]
        add_action(
            expected_relative,
            file.path("posteriors", paste0(new_posterior_name, ".json")),
            pdbtools_rename_json_text(serializable)
        )
        affected[[length(affected) + 1L]] <- list(
            old = old_posterior_name,
            new = new_posterior_name,
            reference = record$reference_posterior_name
        )
    }

    all_references <- vapply(
        posterior_records,
        function(record) {
            reference <- record$reference_posterior_name
            if (is.null(reference) || !length(reference)) "" else as.character(reference[[1L]])
        },
        character(1)
    )
    affected_old <- if (length(affected)) {
        vapply(affected, `[[`, character(1), "old")
    } else {
        character()
    }
    for (item in affected) {
        if (
            identical(item$old, item$new) ||
                is.null(item$reference) ||
                !length(item$reference) ||
                !identical(item$reference, item$new)
        ) {
            next
        }
        users <- which(all_references == item$old)
        if (
            length(users) &&
                any(!vapply(
                    posterior_records[users],
                    function(record) record$name %in% affected_old,
                    logical(1)
                ))
        ) {
            stop(
                "Reference posterior '",
                item$old,
                "' is shared with an unaffected posterior.",
                call. = FALSE
            )
        }
        actions <- pdbtools_rename_add_reference_actions(
            actions,
            root,
            item$old,
            item$new
        )
    }

    alias_path <- file.path(root, "alias", "posteriors.json")
    if (file.exists(alias_path) && length(affected)) {
        aliases <- pdbtools_rename_read_json(alias_path)
        changed <- FALSE
        for (item in affected) {
            for (i in seq_along(aliases)) {
                if (identical(aliases[[i]], item$old)) {
                    aliases[[i]] <- item$new
                    changed <- TRUE
                }
            }
        }
        if (changed) {
            add_action(
                file.path("alias", "posteriors.json"),
                file.path("alias", "posteriors.json"),
                pdbtools_rename_json_text(aliases)
            )
        }
    }

    pdbtools_rename_commit(actions, root)
    list(
        type = type,
        old_name = old_name,
        new_name = new_name,
        posterior_renames = affected,
        affected_posteriors = if (length(affected)) {
            vapply(affected, `[[`, character(1), "new")
        } else {
            character()
        },
        files = vapply(unname(actions), `[[`, character(1), "target")
    )
}

PDBEntryBuilder <- R6::R6Class(
    "PDBEntryBuilder",

    public = list(
        path = NULL,
        adder = NULL,
        data = NULL,
        data_name = NULL,
        stan_file = NULL,
        model_code = NULL,
        model_name = NULL,
        posterior = NULL,
        rp = NULL,
        initialize = function(
            path,
            added_by = "Gerald Press",
            auto_write = TRUE,
            n_threads = 2,
            detect_cores = TRUE,
            n_cores = NULL
        ) {
            required_packages <- c(
                "rstan",
                "posterior",
                "posteriordb",
                "checkmate"
            )
            missing_packages <- required_packages[
                !vapply(
                    required_packages,
                    requireNamespace,
                    logical(1),
                    quietly = TRUE
                )
            ]
            if (length(missing_packages) > 0L) {
                stop(
                    "Missing required packages: ",
                    paste(missing_packages, collapse = ", "),
                    call. = FALSE
                )
            }
            rstan::rstan_options(auto_write = auto_write)
            if (detect_cores) {
                options(mc.cores = parallel::detectCores())
            } else if (!is.null(n_cores)) {
                options(mc.cores = n_cores)
            }
            Sys.setenv(STAN_NUM_THREADS = n_threads)
            self$path <- normalizePath(path, mustWork = TRUE)
            self$adder <- added_by
            private$pdb <- posteriordb::pdb_local(path = self$path)
            invisible(self)
        },
        rename_pdb = function(old_name, new_name, type = NULL) {
            migration <- pdbtools_rename_entity(
                database_path = self$path,
                old_name = old_name,
                new_name = new_name,
                type = type
            )

            # Keep the builder's in-memory state usable after the migration.
            if (identical(self$data_name, old_name)) {
                self$data_name <- new_name
            }
            if (identical(self$model_name, old_name)) {
                self$model_name <- new_name
            }
            if (!is.null(self$stan_file) &&
                identical(basename(self$stan_file), paste0(old_name, ".stan"))) {
                self$stan_file <- file.path(
                    dirname(self$stan_file),
                    paste0(new_name, ".stan")
                )
            }
            if (is.list(self$posterior)) {
                if (identical(self$posterior$data_name, old_name)) {
                    self$posterior$data_name <- new_name
                }
                if (identical(self$posterior$model_name, old_name)) {
                    self$posterior$model_name <- new_name
                }
                for (rename in migration$posterior_renames) {
                    if (identical(self$posterior$name, rename$old)) {
                        self$posterior$name <- rename$new
                    }
                    if (identical(
                        self$posterior$reference_posterior_name,
                        rename$old
                    )) {
                        self$posterior$reference_posterior_name <- rename$new
                    }
                }
            }
            self$refresh()
            invisible(migration)
        },
        link_reference_posterior = function(
            posterior_name,
            reference_posterior_name = posterior_name,
            verify = TRUE
        ) {
            checkmate::assert_string(posterior_name, min.chars = 1L)
            checkmate::assert_string(
                reference_posterior_name,
                min.chars = 1L
            )
            checkmate::assert_flag(verify)

            info_path <- self$get_rpi_path(reference_posterior_name)
            draws_path <- self$get_rp_path(reference_posterior_name)
            if (!file.exists(info_path) || !file.exists(draws_path)) {
                stop(
                    "Cannot link reference posterior `",
                    reference_posterior_name,
                    "`: both its info JSON and draw archive must exist.",
                    call. = FALSE
                )
            }
            reference_info <- jsonlite::read_json(
                info_path,
                simplifyVector = FALSE
            )
            if (!identical(reference_info$name, reference_posterior_name)) {
                stop(
                    "Reference-posterior info name does not match its file name.",
                    call. = FALSE
                )
            }

            posterior_object <- posteriordb::posterior(
                posterior_name,
                private$pdb
            )
            current_reference <- posterior_object$reference_posterior_name
            if (
                !is.null(current_reference) &&
                    !identical(current_reference, reference_posterior_name)
            ) {
                stop(
                    "Posterior already points to a different reference posterior: ",
                    current_reference,
                    call. = FALSE
                )
            }

            changed <- is.null(current_reference)
            if (changed) {
                posterior_object$reference_posterior_name <-
                    reference_posterior_name
                # This is an intentional metadata update to an existing
                # posterior, not replacement of the model/data definition.
                posteriordb::write_pdb(
                    posterior_object,
                    private$pdb,
                    overwrite = TRUE
                )
                if (exists(
                    "pdb_clear_cache",
                    envir = asNamespace("posteriordb"),
                    inherits = FALSE
                )) {
                    try(
                        posteriordb:::pdb_clear_cache(private$pdb),
                        silent = TRUE
                    )
                }
                self$refresh()
            }

            if (verify) {
                posterior_path <- file.path(
                    self$get_pdb_path(),
                    "posterior_database",
                    "posteriors",
                    paste0(posterior_name, ".json")
                )
                written <- jsonlite::read_json(
                    posterior_path,
                    simplifyVector = FALSE
                )
                if (!identical(
                    written$reference_posterior_name,
                    reference_posterior_name
                )) {
                    stop(
                        "Posterior reference link failed round-trip verification.",
                        call. = FALSE
                    )
                }
            }
            invisible(changed)
        },
        link_reference_posterior_from_stan_fit = function(
            stan_fit,
            verify = TRUE
        ) {
            reference_info <- self$get_reference_info(stan_fit)
            if (
                is.null(reference_info) ||
                    is.null(reference_info$name) ||
                    is.null(stan_fit@model_name)
            ) {
                stop(
                    "The Stan fit must contain a model name and reference-posterior info.",
                    call. = FALSE
                )
            }
            self$link_reference_posterior(
                posterior_name = stan_fit@model_name,
                reference_posterior_name = reference_info$name,
                verify = verify
            )
        },
        create_data = function(data, info) {
            if (!inherits(info, "pdb_data_info")) {
                if (!is.list(info)) {
                    stop(
                        "info must be a list or a `pdb_data_info` object",
                        call. = FALSE
                    )
                }
                info <- posteriordb::as.pdb_data_info(info)
            }
            if (!is.list(data)) {
                stop(
                    "data must be a named list containing Stan data",
                    call. = FALSE
                )
            }
            data_names <- names(data)
            if (
                any(is.null(data_names)) ||
                    any(is.na(data_names)) ||
                    any(data_names == "")
            ) {
                stop(
                    "`data` must be a fully named list.",
                    call. = FALSE
                )
            }
            if (anyDuplicated(data_names)) {
                stop(
                    "`data` must not contain duplicate names.",
                    call. = FALSE
                )
            }
            self$data <- data
            posteriordb::as.pdb_data(data, info = info)
        },
        add_data = function(data, info, overwrite = FALSE) {
            pdb_data <- self$create_data(data, info)
            posteriordb::write_pdb(pdb_data, private$pdb, overwrite = overwrite)
            self$refresh()
            invisible(pdb_data)
        },
        generate_model_dims = function(
            model_code,
            data = self$data,
            include = NULL,
            exclude = NULL
        ) {
            always_exclude <- "lp__"

            if (is.null(data)) {
                stop(
                    paste(
                        "Stan data are not stored in this object",
                        "and were not supplied to the method."
                    ),
                    call. = FALSE
                )
            }

            if (!is.list(data)) {
                stop(
                    "`data` must be a named list or `pdb_data` object.",
                    call. = FALSE
                )
            }

            conflicts <- intersect(include, exclude)

            if (length(conflicts) > 0L) {
                stop(
                    paste0(
                        "Variables cannot be both included and excluded: ",
                        paste(conflicts, collapse = ", ")
                    ),
                    call. = FALSE
                )
            }

            forbidden <- intersect(include, always_exclude)

            if (length(forbidden) > 0L) {
                stop(
                    paste0(
                        "These variables cannot be included: ",
                        paste(forbidden, collapse = ", ")
                    ),
                    call. = FALSE
                )
            }

            fit <- suppressWarnings(
                rstan::stan(
                    model_code = as.character(model_code),
                    chains = 1,
                    iter = 4,
                    warmup = 2,
                    refresh = 0,
                    data = data
                )
            )

            model_dims <- fit@par_dims
            available <- names(model_dims)

            if (!is.null(include)) {
                missing_include <- setdiff(
                    include,
                    available
                )

                if (length(missing_include) > 0L) {
                    stop(
                        paste0(
                            "Variables requested for inclusion were not found: ",
                            paste(missing_include, collapse = ", ")
                        ),
                        call. = FALSE
                    )
                }
            }

            if (!is.null(exclude)) {
                missing_exclude <- setdiff(
                    exclude,
                    available
                )

                if (length(missing_exclude) > 0L) {
                    stop(
                        paste0(
                            "Variables requested for exclusion were not found: ",
                            paste(missing_exclude, collapse = ", ")
                        ),
                        call. = FALSE
                    )
                }
            }

            if (is.null(include)) {
                selected <- available
            } else {
                selected <- include
            }

            selected <- setdiff(
                selected,
                c(exclude, always_exclude)
            )

            if (length(selected) == 0L) {
                stop(
                    "Parameter selection produced no variables.",
                    call. = FALSE
                )
            }

            model_dims <- model_dims[selected]

            invisible(model_dims)
        },
        create_model_code = function(stan_file, info) {
            if (!inherits(info, "pdb_model_info")) {
                if (!is.list(info)) {
                    stop(
                        "info must be a list or a `pdb_model_info` object",
                        call. = FALSE
                    )
                }
                info <- posteriordb::as.pdb_model_info(info)
            }
            stan_model <- tryCatch(
                {
                    rstan::stan_model(file = stan_file)
                },
                error = function(e) {
                    stop(
                        "Stan model could not be compiled:",
                        conditionMessage(e),
                        call. = FALSE
                    )
                }
            )
            self$stan_file <- stan_file
            pdb_model_code <- posteriordb::as.pdb_model_code(
                stan_model,
                info = info
            )
            self$model_code <- pdb_model_code
            pdb_model_code
        },
        add_model_code = function(stan_file, info, overwrite = FALSE) {
            # Uses private$pdb
            pdb_model_code <- self$create_model_code(stan_file, info)
            posteriordb::write_pdb(
                pdb_model_code,
                private$pdb,
                overwrite = overwrite
            )
            self$refresh()
            invisible(pdb_model_code)
        },
        create_posterior = function(
            info,
            include = NULL,
            exclude = c("lp__")
        ) {
            # Uses private$pdb
            metadata_names <- intersect(
                c("keywords", "urls", "references"),
                names(info)
            )
            metadata <- info[metadata_names]

            if (!inherits(info, "pdb_posterior")) {
                if (!is.list(info)) {
                    stop(
                        "Posterior must be a list or `pdb_posterior` object.",
                        call. = FALSE
                    )
                }
                if (is.null(info$dimensions)) {
                    model_code <- info$pdb_model_code
                    if (is.null(model_code)) {
                        model_code <- self$model_code
                    }
                    if (is.null(model_code)) {
                        model_code <- self$stan_file
                    }

                    data <- info$pdb_data
                    if (is.null(data)) {
                        data <- self$data
                    }

                    if (is.null(model_code)) {
                        stop(
                            paste(
                                "Posterior dimensions are missing and no Stan",
                                "model code is available to generate them."
                            ),
                            call. = FALSE
                        )
                    }
                    message("Dimensions not provided. Attemping to generate...")
                    info[["dimensions"]] <- self$generate_model_dims(
                        model_code = model_code,
                        data = data,
                        include = include,
                        exclude = exclude
                    )
                    #Testing: product of dimensions for reference computation
                    # if (any(sapply(info[["dimensions"]], length) > 1)) {
                    #     message(
                    #         "Converting multi-dimensional parameters to one dimension..."
                    #     )
                    #     info[["dimensions"]] <- lapply(
                    #         info[["dimensions"]],
                    #         prod
                    #     )
                    # }
                }
                info <- posteriordb::as.pdb_posterior(
                    info,
                    pdb = private$pdb
                )
            } else {
                # Rebuild the object so it is associated with this local PDB.
                info <- posteriordb::as.pdb_posterior(
                    unclass(info),
                    pdb = private$pdb
                )
            }

            # posteriordb 0.3.6 drops optional posterior metadata while
            # coercing; add those fields back before writing.
            for (field in names(metadata)) {
                if (!is.null(metadata[[field]])) {
                    info[[field]] <- metadata[[field]]
                }
            }

            invisible(info)
        },
        prepare_posterior = function(
            data_name,
            model_name,
            keywords = NULL,
            references = NULL,
            framework = "stan",
            dimensions = NULL,
            added_date = Sys.Date(),
            added_by = "Gerald Press",
            parameters_include = NULL,
            parameters_exclude = NULL
        ) {
            # Uses private$pdb
            pdb_data <- posteriordb::pdb_data(data_name, private$pdb)
            pdb_model_info <- posteriordb::model_info(model_name, private$pdb)
            pdb_model_code <- posteriordb::model_code(
                model_name,
                framework = framework,
                pdb = private$pdb
            )
            if (!inherits(pdb_data, "pdb_data")) {
                stop("Data must be a `pdb_data` object.", call. = FALSE)
            }
            if (!inherits(pdb_model_info, "pdb_model_info")) {
                stop(
                    "Model info must be a `pdb_model_info` object.",
                    call. = FALSE
                )
            }
            if (!inherits(pdb_model_code, "pdb_model_code")) {
                stop(
                    "Model code must be a `pdb_model_code` object.",
                    call. = FALSE
                )
            }
            posterior <- list(
                pdb_model_code = pdb_model_code,
                pdb_data = pdb_data,
                keywords = keywords,
                references = references,
                framework = framework,
                dimensions = dimensions,
                added_date = added_date,
                added_by = added_by
            )
            if (is.null(posterior$dimensions)) {
                message("Dimensions not provided. Attemping to generate...")
                posterior$dimensions <- self$generate_model_dims(
                    model_code = pdb_model_code,
                    data = pdb_data,
                    include = parameters_include,
                    exclude = parameters_exclude
                )
                # if (any(sapply(posterior$dimensions, length) > 1)) {
                #     message(
                #         "Converting multi-dimensional parameters to one dimension..."
                #     )
                #     posterior$dimensions <- lapply(
                #         posterior$dimensions,
                #         prod
                #     )
                # }
            }
            po <- self$create_posterior(posterior)
            self$set_posterior(po)
            invisible(po)
        },
        check_posterior = function(p = self$posterior, checks = FALSE) {
            if (is.null(p)) {
                stop(
                    "Posterior not stored in object, or defined in function call.",
                    call. = FALSE
                )
            }
            if (!inherits(p, "pdb_posterior")) {
                stop(
                    "Posterior must be a `pdb_posterior` object.",
                    call. = FALSE
                )
            }
            if (
                !is.logical(checks) ||
                    length(checks) != 1L ||
                    is.na(checks)
            ) {
                stop("`checks` must be TRUE or FALSE.", call. = FALSE)
            }
            posteriordb::check_pdb_posterior(
                p,
                run_stan_code_checks = checks
            )
            invisible(p)
        },
        add_posterior = function(
            spec,
            overwrite = FALSE,
            dry_run = TRUE
        ) {
            if (
                !is.logical(overwrite) ||
                    length(overwrite) != 1L ||
                    is.na(overwrite)
            ) {
                stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
            }
            if (
                !is.logical(dry_run) ||
                    length(dry_run) != 1L ||
                    is.na(dry_run)
            ) {
                stop("`dry_run` must be TRUE or FALSE.", call. = FALSE)
            }

            if (inherits(spec, "pdb_posterior")) {
                po <- self$create_posterior(spec)
            } else {
                if (!is.list(spec)) {
                    stop(
                        paste(
                            "`spec` must be a list or a",
                            "`pdb_posterior` object."
                        ),
                        call. = FALSE
                    )
                }

                spec_names <- names(spec)
                if (
                    is.null(spec_names) ||
                        any(is.na(spec_names)) ||
                        any(spec_names == "")
                ) {
                    stop("`spec` must be a fully named list.", call. = FALSE)
                }
                if (anyDuplicated(spec_names)) {
                    stop(
                        "`spec` must not contain duplicate names.",
                        call. = FALSE
                    )
                }

                uses_database_names <- all(
                    c("data_name", "model_name") %in% spec_names
                )

                if (uses_database_names) {
                    allowed <- names(formals(self$prepare_posterior))
                    unknown <- setdiff(spec_names, allowed)
                    if (length(unknown) > 0L) {
                        stop(
                            paste0(
                                "Unknown posterior specification field(s): ",
                                paste(unknown, collapse = ", ")
                            ),
                            call. = FALSE
                        )
                    }
                    po <- do.call(self$prepare_posterior, spec)
                } else {
                    required <- c("pdb_model_code", "pdb_data")
                    missing <- setdiff(required, spec_names)
                    if (length(missing) > 0L) {
                        stop(
                            paste0(
                                "`spec` must contain either `data_name` and ",
                                "`model_name`, or: ",
                                paste(required, collapse = ", "),
                                ". Missing: ",
                                paste(missing, collapse = ", ")
                            ),
                            call. = FALSE
                        )
                    }

                    include <- spec$parameters_include
                    exclude <- spec$parameters_exclude
                    spec[c(
                        "parameters_include",
                        "parameters_exclude"
                    )] <- NULL

                    po <- self$create_posterior(
                        info = spec,
                        include = include,
                        exclude = exclude
                    )
                }
            }

            self$set_posterior(po)

            if (dry_run) {
                return(invisible(po))
            }

            posteriordb::write_pdb(
                po,
                private$pdb,
                overwrite = overwrite
            )
            self$refresh()

            po <- posteriordb::posterior(po$name, private$pdb)
            self$check_posterior(po, checks = TRUE)
            self$set_posterior(po)
            self$refresh()
            invisible(po)
        },
        run_entry_workflow = function(
            entry_spec,
            model = NULL,
            register = TRUE,
            sample = TRUE,
            write = FALSE,
            overwrite = FALSE
        ) {
            results <- self$run_workflows(
                model = model,
                entries = list(entry = entry_spec),
                register = register,
                sample = sample,
                write = write,
                overwrite = overwrite,
                continue_on_error = FALSE
            )
            invisible(results[[1L]])
        },
        run_workflows = function(
            model = NULL,
            entries,
            register = TRUE,
            sample = TRUE,
            write = FALSE,
            overwrite = FALSE,
            continue_on_error = TRUE,
            save_failed_fits = FALSE,
            failed_fit_dir = NULL
        ) {
            flags <- list(
                register = register,
                sample = sample,
                write = write,
                overwrite = overwrite,
                continue_on_error = continue_on_error,
                save_failed_fits = save_failed_fits
            )
            for (flag_name in names(flags)) {
                checkmate::assert_flag(flags[[flag_name]])
            }
            if (write && !sample) {
                stop("`write = TRUE` requires `sample = TRUE`.", call. = FALSE)
            }
            if (save_failed_fits) {
                if (
                    is.null(failed_fit_dir) ||
                        !is.character(failed_fit_dir) ||
                        length(failed_fit_dir) != 1L ||
                        is.na(failed_fit_dir) ||
                        !nzchar(failed_fit_dir)
                ) {
                    stop(
                        "`failed_fit_dir` must be a non-empty path when `save_failed_fits = TRUE`.",
                        call. = FALSE
                    )
                }
                dir.create(
                    failed_fit_dir,
                    recursive = TRUE,
                    showWarnings = FALSE
                )
                if (!dir.exists(failed_fit_dir)) {
                    stop(
                        "Could not create `failed_fit_dir`.",
                        call. = FALSE
                    )
                }
            }

            save_failed_fit <- function(fit, posterior_name, reason) {
                if (!save_failed_fits || !inherits(fit, "stanfit")) {
                    return(NULL)
                }

                safe_name <- gsub(
                    "[^A-Za-z0-9_.-]",
                    "_",
                    posterior_name
                )
                stamp <- format(
                    Sys.time(),
                    "%Y%m%d_%H%M%S"
                )
                prefix <- file.path(
                    failed_fit_dir,
                    paste0(
                        safe_name,
                        "_",
                        stamp,
                        "_pid",
                        Sys.getpid()
                    )
                )
                fit_path <- paste0(prefix, ".rds")
                diagnostic_path <- paste0(prefix, ".diagnostics.rds")

                saved <- tryCatch(
                    {
                        saveRDS(fit, fit_path, compress = "gzip")
                        saveRDS(
                            list(
                                posterior_name = posterior_name,
                                saved_at = Sys.time(),
                                reason = reason,
                                reference_info = self$get_reference_info(fit)
                            ),
                            diagnostic_path,
                            compress = "gzip"
                        )
                        message(
                            "Saved failed fit for ",
                            posterior_name,
                            " to ",
                            fit_path
                        )
                        fit_path
                    },
                    error = function(e) {
                        warning(
                            paste0(
                                "Could not save failed fit for ",
                                posterior_name,
                                ": ",
                                conditionMessage(e)
                            ),
                            call. = FALSE
                        )
                        NULL
                    }
                )
                saved
            }
            if (!is.list(entries) || length(entries) == 0L) {
                stop("`entries` must be a non-empty list.", call. = FALSE)
            }
            if (is.null(names(entries))) {
                names(entries) <- paste0("entry_", seq_along(entries))
            }
            if (anyNA(names(entries)) || any(names(entries) == "")) {
                stop(
                    "Every entry specification must have a name.",
                    call. = FALSE
                )
            }
            if (anyDuplicated(names(entries))) {
                stop("Entry specification names must be unique.", call. = FALSE)
            }

            common_model_name <- NULL
            model_object <- NULL
            if (!is.null(model)) {
                if (!is.list(model)) {
                    stop("`model` must be NULL or a list.", call. = FALSE)
                }
                missing_model_fields <- setdiff(
                    c("stan_file", "info"),
                    names(model)
                )
                if (length(missing_model_fields) > 0L) {
                    stop(
                        paste0(
                            "The model specification is missing: ",
                            paste(missing_model_fields, collapse = ", ")
                        ),
                        call. = FALSE
                    )
                }
                common_model_name <- model$info$name
                if (is.null(common_model_name)) {
                    stop("`model$info$name` is required.", call. = FALSE)
                }
                if (register) {
                    model_overwrite <- if (is.null(model$overwrite)) {
                        overwrite
                    } else {
                        model$overwrite
                    }
                    checkmate::assert_flag(model_overwrite)
                    model_object <- self$add_model_code(
                        stan_file = model$stan_file,
                        info = model$info,
                        overwrite = model_overwrite
                    )
                }
            }

            results <- vector("list", length(entries))
            names(results) <- names(entries)

            for (i in seq_along(entries)) {
                entry_name <- names(entries)[[i]]
                spec <- entries[[i]]
                fit <- NULL
                posterior_object <- NULL
                posterior_name <- entry_name

                message(
                    "Starting workflow ",
                    i,
                    " of ",
                    length(entries),
                    ": ",
                    entry_name
                )

                result <- tryCatch(
                    {
                        if (!is.list(spec)) {
                            stop("The entry specification must be a list.")
                        }
                        data_spec <- spec$data
                        posterior_spec <- spec$posterior
                        reference_spec <- spec$reference

                        if (is.null(posterior_spec)) {
                            posterior_spec <- list()
                        }
                        if (!is.list(posterior_spec)) {
                            stop("`entry_spec$posterior` must be a list.")
                        }

                        data_name <- posterior_spec$data_name
                        if (is.null(data_name) && !is.null(data_spec)) {
                            data_name <- data_spec$info$name
                        }
                        model_name <- posterior_spec$model_name
                        if (is.null(model_name)) {
                            model_name <- common_model_name
                        }
                        if (is.null(data_name) || is.null(model_name)) {
                            stop(
                                "Each entry must resolve `data_name` and `model_name`."
                            )
                        }

                        entry_overwrite <- if (is.null(spec$overwrite)) {
                            overwrite
                        } else {
                            spec$overwrite
                        }
                        checkmate::assert_flag(entry_overwrite)

                        if (register) {
                            if (
                                !is.list(data_spec) ||
                                    is.null(data_spec$data) ||
                                    is.null(data_spec$info)
                            ) {
                                stop(
                                    paste0(
                                        "`entry_spec$data` must contain `data` ",
                                        "and `info` when `register = TRUE`."
                                    )
                                )
                            }
                            self$add_data(
                                data = data_spec$data,
                                info = data_spec$info,
                                overwrite = entry_overwrite
                            )

                            posterior_spec$data_name <- data_name
                            posterior_spec$model_name <- model_name
                            allowed_posterior_args <- names(
                                formals(self$prepare_posterior)
                            )
                            unknown_posterior_args <- setdiff(
                                names(posterior_spec),
                                allowed_posterior_args
                            )
                            if (length(unknown_posterior_args) > 0L) {
                                stop(
                                    paste0(
                                        "Unknown posterior specification fields: ",
                                        paste(
                                            unknown_posterior_args,
                                            collapse = ", "
                                        )
                                    )
                                )
                            }
                            posterior_object <- do.call(
                                self$prepare_posterior,
                                posterior_spec
                            )
                            posterior_object <- self$add_posterior(
                                spec = posterior_object,
                                overwrite = entry_overwrite,
                                dry_run = FALSE
                            )
                        }

                        posterior_name <- if (!is.null(spec$posterior_name)) {
                            spec$posterior_name
                        } else if (!is.null(posterior_object$name)) {
                            posterior_object$name
                        } else {
                            paste(data_name, model_name, sep = "-")
                        }

                        checks <- NULL
                        failed_checks <- character()
                        failed_required_checks <- character()
                        divergences_by_chain <- NULL
                        total_divergences <- NA_real_
                        written <- FALSE
                        info_path <- NULL
                        draws_path <- NULL
                        summary_paths <- NULL
                        failed_fit_path <- NULL

                        if (sample) {
                            if (
                                !is.list(reference_spec) ||
                                    is.null(reference_spec$sampling_args)
                            ) {
                                stop(
                                    paste0(
                                        "`entry_spec$reference$sampling_args` ",
                                        "is required when `sample = TRUE`."
                                    )
                                )
                            }
                            reference_args <- reference_spec
                            reference_args$sampling_args <- NULL
                            allowed_reference_args <- setdiff(
                                names(formals(self$compute_reference_draws)),
                                c(
                                    "posterior_name",
                                    "sampling_args",
                                    "auto_check",
                                    "write",
                                    "overwrite"
                                )
                            )
                            unknown_reference_args <- setdiff(
                                names(reference_args),
                                allowed_reference_args
                            )
                            if (length(unknown_reference_args) > 0L) {
                                stop(
                                    paste0(
                                        "Unknown reference specification fields: ",
                                        paste(
                                            unknown_reference_args,
                                            collapse = ", "
                                        )
                                    )
                                )
                            }
                            compute_args <- c(
                                list(
                                    posterior_name = posterior_name,
                                    sampling_args = reference_spec$sampling_args,
                                    auto_check = FALSE,
                                    write = FALSE,
                                    overwrite = entry_overwrite
                                ),
                                reference_args
                            )
                            fit <- do.call(
                                self$compute_reference_draws,
                                compute_args
                            )

                            fit_info <- self$get_reference_info(fit)
                            divergences_by_chain <-
                                fit_info$diagnostics$divergent_transitions
                            total_divergences <- if (
                                length(divergences_by_chain) > 0L &&
                                    !anyNA(divergences_by_chain) &&
                                    all(is.finite(divergences_by_chain))
                            ) {
                                sum(divergences_by_chain)
                            } else {
                                NA_real_
                            }

                            if (
                                is.na(total_divergences) ||
                                    total_divergences > 0L
                            ) {
                                failed_checks <- "divergent_transitions"
                                failed_required_checks <- failed_checks
                            } else {
                                checks <- self$get_checks_from_stanfit(fit)
                                failed_checks <- names(checks)[
                                    !vapply(checks, isTRUE, logical(1))
                                ]
                                required_checks <- c(
                                    "ndraws_is_10k",
                                    "nchains_is_gte_4",
                                    "r_hat_below_1_01",
                                    "efmi_above_0_2",
                                    "abs_mean_lag1_ac_below_0_05"
                                )
                                failed_required_checks <- intersect(
                                    failed_checks,
                                    required_checks
                                )
                            }

                            if (length(failed_required_checks) > 0L) {
                                failed_fit_path <- save_failed_fit(
                                    fit,
                                    posterior_name,
                                    list(
                                        status = "failed_checks",
                                        failed_checks = failed_checks,
                                        failed_required_checks = failed_required_checks,
                                        total_divergences = total_divergences
                                    )
                                )
                            }

                            if (length(failed_required_checks) == 0L) {
                                fit <- self$check_draws_from_stanfit(fit)
                                if (write) {
                                    info_path <- self$write_rpi_from_stan_fit(
                                        fit,
                                        overwrite = entry_overwrite,
                                        verify = TRUE
                                    )
                                    draws_path <- self$write_rpd_from_stan_fit(
                                        fit,
                                        overwrite = entry_overwrite,
                                        verify = TRUE
                                    )
                                    summary_paths <- self$write_summary_statistics_from_stan_fit(
                                        fit,
                                        overwrite = entry_overwrite,
                                        verify = TRUE
                                    )
                                    self$verify_reference_files(fit)
                                    self$link_reference_posterior_from_stan_fit(fit)
                                    written <- TRUE
                                }
                            }
                        }

                        status <- if (length(failed_required_checks) > 0L) {
                            "failed_checks"
                        } else {
                            "completed"
                        }
                        message(entry_name, " finished with status: ", status)
                        list(
                            name = entry_name,
                            posterior_name = posterior_name,
                            status = status,
                            model = model_object,
                            posterior = posterior_object,
                            fit = fit,
                            diagnostics = if (is.null(fit)) {
                                NULL
                            } else {
                                self$get_reference_info(fit)$diagnostics
                            },
                            divergences_by_chain = divergences_by_chain,
                            total_divergences = total_divergences,
                            checks = checks,
                            failed_checks = failed_checks,
                            failed_required_checks = failed_required_checks,
                            written = written,
                            info_path = info_path,
                            draws_path = draws_path,
                            summary_paths = summary_paths,
                            failed_fit_path = failed_fit_path,
                            error = NULL
                        )
                    },
                    error = function(e) {
                        failed_fit_path <- save_failed_fit(
                            fit,
                            if (is.null(posterior_name)) entry_name else posterior_name,
                            list(
                                status = "error",
                                error = conditionMessage(e)
                            )
                        )
                        message(entry_name, " failed: ", conditionMessage(e))
                        list(
                            name = entry_name,
                            posterior_name = NULL,
                            status = "error",
                            model = model_object,
                            posterior = posterior_object,
                            fit = fit,
                            diagnostics = NULL,
                            divergences_by_chain = NULL,
                            total_divergences = NA_real_,
                            checks = NULL,
                            failed_checks = character(),
                            failed_required_checks = character(),
                            written = FALSE,
                            info_path = NULL,
                            draws_path = NULL,
                            summary_paths = NULL,
                            failed_fit_path = failed_fit_path,
                            error = conditionMessage(e)
                        )
                    }
                )
                results[[i]] <- result

                if (
                    identical(result$status, "error") &&
                        !continue_on_error
                ) {
                    stop(result$error, call. = FALSE)
                }
            }

            class(results) <- c("pdb_workflow_results", "list")
            invisible(results)
        },
        summarize_workflow_results = function(results) {
            if (!inherits(results, "pdb_workflow_results")) {
                stop(
                    "`results` must come from `run_workflows()`.",
                    call. = FALSE
                )
            }
            data.frame(
                name = names(results),
                posterior = vapply(
                    results,
                    function(x) {
                        if (is.null(x$posterior_name)) {
                            NA_character_
                        } else {
                            x$posterior_name
                        }
                    },
                    character(1)
                ),
                status = vapply(results, function(x) x$status, character(1)),
                divergences = vapply(
                    results,
                    function(x) x$total_divergences,
                    numeric(1)
                ),
                failed_checks = vapply(
                    results,
                    function(x) paste(x$failed_checks, collapse = ", "),
                    character(1)
                ),
                failed_required_checks = vapply(
                    results,
                    function(x) {
                        paste(x$failed_required_checks, collapse = ", ")
                    },
                    character(1)
                ),
                written = vapply(
                    results,
                    function(x) x$written,
                    logical(1)
                ),
                error = vapply(
                    results,
                    function(x) {
                        if (is.null(x$error)) NA_character_ else x$error
                    },
                    character(1)
                ),
                row.names = NULL
            )
        },
        compute_reference_draws = function(
            posterior_name,
            sampling_args,
            control_args = NULL, #must be a list
            comments = NULL,
            added_by = self$adder,
            auto_check = TRUE,
            write = FALSE,
            overwrite = FALSE
        ) {
            if (!auto_check && write) {
                stop("Can't write draws without checking them first.")
            }
            posterior_objects <- self$get_posterior_modeldata(posterior_name)
            sa <- list(
                object = posterior_objects$stan_model,
                data = posterior_objects$data
            )
            if (is.null(control_args)) {
                sa <- c(sa, sampling_args)
                inference <- sampling_args
            } else if (inherits(control_args, "list")) {
                sa <- c(sa, sampling_args, control_args)
                inference <- c(sampling_args, control_args)
            } else {
                stop("`control_args` must be null or a list.", call. = FALSE)
            }
            stan_fit <- do.call(rstan::sampling, sa)
            stan_fit@model_name <- posterior_name
            stan_fit <- self$set_reference_info(
                stan_fit,
                list(
                    name = posterior_name,
                    inference = list(
                        method = "stan_sampling",
                        method_arguments = inference
                    ),
                    diagnostics = self$get_diagnostics(stan_fit),
                    checks_made = list(),
                    comments = comments,
                    added_by = added_by,
                    added_date = Sys.Date(),
                    versions = private$get_sampling_version_info()
                )
            )
            if (auto_check) {
                stan_fit <- self$check_draws_from_stanfit(stan_fit)
            }
            if (write) {
                self$write_rpi_from_stan_fit(stan_fit, overwrite = overwrite)
                self$write_rpd_from_stan_fit(
                    stan_fit,
                    overwrite = overwrite,
                    verify = TRUE
                )
                self$write_summary_statistics_from_stan_fit(
                    stan_fit,
                    overwrite = overwrite,
                    verify = TRUE
                )
                self$verify_reference_files(stan_fit)
                self$link_reference_posterior_from_stan_fit(stan_fit)
            }
            invisible(stan_fit)
        },
        write_rpi_from_stan_fit = function(
            stan_fit,
            overwrite = FALSE,
            verify = TRUE
        ) {
            if (is.null(self$get_reference_info(stan_fit))) {
                stop(
                    "`stan_fit` object must contain reference draw info. Call `add_rpi_from_stanfit()`",
                    call. = FALSE
                )
            }
            info <- self$get_reference_info(stan_fit)
            if (is.null(info$checks_made)) {
                stop(
                    "Draws must be checked before writing. Call `check_draws_from_stanfit()`",
                    call. = FALSE
                )
            }
            required_checks <- c(
                "ndraws_is_10k",
                "nchains_is_gte_4",
                "r_hat_below_1_01",
                "efmi_above_0_2",
                "abs_mean_lag1_ac_below_0_05"
            )
            failed_checks <- required_checks[
                !vapply(
                    info$checks_made[required_checks],
                    isTRUE,
                    logical(1)
                )
            ]
            if (length(failed_checks) > 0L) {
                stop(
                    paste0(
                        "The following checks did not pass: ",
                        paste(failed_checks, collapse = ", ")
                    ),
                    call. = FALSE
                )
            }
            info_path <- self$get_rpi_path(stan_fit@model_name)
            dir.create(
                dirname(info_path),
                recursive = TRUE,
                showWarnings = FALSE
            )
            if (file.exists(info_path) && !overwrite) {
                stop(
                    "Reference-posterior information already exists.",
                    call. = FALSE
                )
            }
            temp_info <- tempfile(
                pattern = paste0(stan_fit@model_name, "-"),
                tmpdir = dirname(info_path),
                fileext = ".json"
            )
            on.exit(unlink(temp_info, force = TRUE), add = TRUE)
            jsonlite::write_json(
                info,
                temp_info,
                pretty = TRUE,
                auto_unbox = TRUE,
                null = "null",
                digits = NA
            )
            if (!file.copy(temp_info, info_path, overwrite = overwrite)) {
                stop(
                    "Failed to write reference-posterior information.",
                    call. = FALSE
                )
            }
            if (verify) {
                loaded_info <- jsonlite::read_json(
                    info_path,
                    simplifyVector = FALSE
                )
                if (!identical(loaded_info$name, info$name)) {
                    stop(
                        "Written reference-posterior information failed verification.",
                        call. = FALSE
                    )
                }
            }
            invisible(info_path)
        },
        write_rpd_from_stan_fit = function(
            stan_fit,
            overwrite = FALSE,
            verify = TRUE
        ) {
            if (!file.exists(self$get_rpi_path(stan_fit@model_name))) {
                stop(
                    "Make sure to write the reference posterior information to disk before writing the draws",
                    call. = FALSE
                )
            }
            to_keep <- self$get_posterior_dims(stan_fit@model_name)
            draws <- posterior::subset_draws(
                posterior::as_draws_array(stan_fit),
                variable = names(to_keep)
            )
            rp_path <- self$get_rp_path(stan_fit@model_name)
            dir.create(dirname(rp_path), recursive = TRUE, showWarnings = FALSE)
            if (file.exists(rp_path) && !overwrite) {
                stop("Reference-draw archive already exists.", call. = FALSE)
            }
            temp_dir <- tempfile(
                pattern = paste0(stan_fit@model_name, "-"),
                tmpdir = dirname(rp_path)
            )
            dir.create(temp_dir)
            json_path <- file.path(
                temp_dir,
                paste0(stan_fit@model_name, ".json")
            )
            temp_zip <- tempfile(
                pattern = paste0(stan_fit@model_name, "-"),
                tmpdir = dirname(rp_path),
                fileext = ".json.zip"
            )
            on.exit(
                {
                    unlink(temp_zip, force = TRUE)
                    unlink(temp_dir, recursive = TRUE, force = TRUE)
                },
                add = TRUE
            )
            jsonlite::write_json(draws, json_path, digits = NA, null = "null")
            zip_status <- utils::zip(
                zipfile = temp_zip,
                files = json_path,
                flags = "-jq"
            )
            if (!identical(zip_status, 0L) || !file.exists(temp_zip)) {
                stop(
                    "Failed to create the reference-draw archive.",
                    call. = FALSE
                )
            }
            if (!file.copy(temp_zip, rp_path, overwrite = overwrite)) {
                stop(
                    "Failed to write the reference-draw archive.",
                    call. = FALSE
                )
            }
            if (verify) {
                self$verify_reference_files(stan_fit, expected_draws = draws)
            }
            invisible(rp_path)
        },
        read_reference_files = function(posterior_name) {
            info_path <- self$get_rpi_path(posterior_name)
            draws_path <- self$get_rp_path(posterior_name)
            if (!file.exists(info_path) || !file.exists(draws_path)) {
                stop(
                    "Reference-posterior information or draws are missing.",
                    call. = FALSE
                )
            }
            archive <- utils::unzip(draws_path, list = TRUE)
            expected_member <- paste0(posterior_name, ".json")
            if (
                nrow(archive) != 1L ||
                    basename(archive$Name[[1]]) != expected_member
            ) {
                stop(
                    "Reference-draw archive has an unexpected structure.",
                    call. = FALSE
                )
            }
            extraction_dir <- tempfile(pattern = paste0(posterior_name, "-"))
            dir.create(extraction_dir)
            on.exit(
                unlink(extraction_dir, recursive = TRUE, force = TRUE),
                add = TRUE
            )
            extracted <- utils::unzip(
                draws_path,
                files = archive$Name[[1]],
                exdir = extraction_dir
            )
            list(
                info = jsonlite::read_json(info_path, simplifyVector = FALSE),
                draws = jsonlite::read_json(extracted, simplifyVector = TRUE),
                archive_member = archive$Name[[1]]
            )
        },
        verify_reference_files = function(stan_fit, expected_draws = NULL) {
            posterior_name <- stan_fit@model_name
            if (is.null(expected_draws)) {
                variables <- names(self$get_posterior_dims(posterior_name))
                expected_draws <- posterior::subset_draws(
                    posterior::as_draws_array(stan_fit),
                    variable = variables
                )
            }
            loaded <- self$read_reference_files(posterior_name)
            if (!identical(loaded$info$name, posterior_name)) {
                stop(
                    "Reference-posterior name changed during the write/read round trip.",
                    call. = FALSE
                )
            }
            loaded_values <- unlist(
                loaded$draws,
                recursive = TRUE,
                use.names = FALSE
            )
            expected_values <- as.numeric(unclass(expected_draws))
            if (length(loaded_values) != length(expected_values)) {
                stop(
                    "Reference-draw dimensions changed during the write/read round trip.",
                    call. = FALSE
                )
            }
            if (
                !isTRUE(all.equal(
                    as.numeric(loaded_values),
                    expected_values,
                    tolerance = sqrt(.Machine$double.eps),
                    check.attributes = FALSE
                ))
            ) {
                stop(
                    "Reference-draw values changed during the write/read round trip.",
                    call. = FALSE
                )
            }
            invisible(TRUE)
        },
        get_checks_from_stanfit = function(stan_fit, diagnostics = NULL) {
            if (is.null(diagnostics)) {
                diagnostics <- self$get_diagnostics(stan_fit)
            }
            if (
                anyNA(diagnostics$divergent_transitions) ||
                    any(!is.finite(diagnostics$divergent_transitions)) ||
                    any(diagnostics$divergent_transitions > 0)
            ) {
                stop(
                    "There were invalid or divergent transitions during sampling.",
                    call. = FALSE
                )
            }
            ess_failures <- self$get_ess_bounds_failures(diagnostics)
            ess_within_bounds <- ess_failures$total_count == 0L
            diagnostic_draws <- posterior::subset_draws(
                posterior::as_draws_array(stan_fit),
                variable = names(diagnostics$rhat)
            )
            list(
                ndraws_is_10k = diagnostics$ndraws == 10000,
                nchains_is_gte_4 = diagnostics$nchains >= 4,
                r_hat_below_1_01 = self$check_rhat(
                    draws = diagnostic_draws,
                    rhat = diagnostics$rhat
                ),
                ess_within_bounds = ess_within_bounds,
                efmi_above_0_2 = !anyNA(diagnostics$efmi) &&
                    all(is.finite(diagnostics$efmi)) &&
                    all(diagnostics$efmi >= 0.2),
                abs_mean_lag1_ac_below_0_05 = all(is.finite(
                    diagnostics$mean_lag1_ac[
                        !is.na(diagnostics$mean_lag1_ac)
                    ]
                )) &&
                    all(diagnostics$mean_lag1_ac <= 0.05, na.rm = TRUE)
            )
        },
        check_draws_from_stanfit = function(stan_fit) {
            stan_info <- self$get_reference_info(stan_fit)
            if (is.null(stan_info)) {
                stop(
                    "`stan_fit` does not have reference-draw information attached.",
                    call. = FALSE
                )
            }
            if (is.null(stan_info$diagnostics)) {
                stan_info$diagnostics <- self$get_diagnostics(stan_fit)
            }
            checks_made <- self$get_checks_from_stanfit(
                stan_fit,
                stan_info$diagnostics
            )
            required_checks <- c(
                "ndraws_is_10k",
                "nchains_is_gte_4",
                "r_hat_below_1_01",
                "efmi_above_0_2",
                "abs_mean_lag1_ac_below_0_05"
            )
            failed_checks <- required_checks[
                !vapply(checks_made[required_checks], isTRUE, logical(1))
            ]
            if (length(failed_checks) > 0L) {
                stop(
                    paste0(
                        "The following checks did not pass: ",
                        paste(failed_checks, collapse = ", ")
                    ),
                    call. = FALSE
                )
            }
            stan_info$checks_made <- checks_made
            stan_fit <- self$set_reference_info(stan_fit, stan_info)
            invisible(stan_fit)
        },
        generate_ess_bounds = function(ndraws = 10000) {
            approx_ess_sd <- sqrt(7) * sqrt(ndraws)
            bnds <- ndraws + 4 * c(approx_ess_sd, -approx_ess_sd)
            list(ess_bulk = bnds, ess_tail = bnds)
        },
        get_ess_bounds_failures = function(diagnostics) {
            ess_bounds <- self$generate_ess_bounds(diagnostics$ndraws)
            find_failures <- function(x, bounds) {
                failed <- is.na(x) |
                    !is.finite(x) |
                    x < min(bounds) |
                    x > max(bounds)
                names(x)[failed]
            }
            bulk <- find_failures(
                diagnostics$effective_sample_size_bulk,
                ess_bounds$ess_bulk
            )
            tail <- find_failures(
                diagnostics$effective_sample_size_tail,
                ess_bounds$ess_tail
            )
            any <- union(bulk, tail)
            list(
                bulk = bulk,
                tail = tail,
                any = any,
                bulk_count = length(bulk),
                tail_count = length(tail),
                total_count = length(any),
                bounds = ess_bounds
            )
        },
        is_within_bounds = function(x, bounds) {
            length(x) > 0L &&
                !anyNA(x) &&
                all(is.finite(x)) &&
                all(x >= min(bounds) & x <= max(bounds))
        },
        compute_ac = function(x) {
            x <- posterior::as_draws_array(x)
            var_names <- posterior::variables(x)
            abs(vapply(
                var_names,
                function(name) {
                    posterior::autocorrelation(
                        posterior::extract_variable(x, name)
                    )[2]
                },
                numeric(1)
            ))
        },
        compute_mean_lag1_ac = function(x) {
            checkmate::assert_class(x, "draws")

            x <- posterior::as_draws_array(x)

            variable_names <- dimnames(x)$variable
            n_chains <- dim(x)[2]
            n_variables <- dim(x)[3]

            # First reject genuinely invalid values.
            all_finite <- apply(
                x,
                3,
                function(z) all(is.finite(z))
            )

            if (any(!all_finite)) {
                bad_variables <- variable_names[!all_finite]

                stop(
                    paste0(
                        "Non-finite draws were found in: ",
                        paste(head(bad_variables, 20), collapse = ", "),
                        if (length(bad_variables) > 20) " ..." else ""
                    ),
                    call. = FALSE
                )
            }

            # A global constant can reasonably have undefined autocorrelation.
            globally_constant <- apply(
                x,
                3,
                self$is_constant
            )

            # A variable constant in only some chains may indicate a stuck chain.
            constant_by_chain <- apply(
                x,
                c(2, 3),
                self$is_constant
            )

            partially_constant <- apply(
                constant_by_chain,
                2,
                any
            ) &
                !globally_constant

            if (any(partially_constant)) {
                bad_variables <- variable_names[partially_constant]

                stop(
                    paste0(
                        "Some variables were constant in only a subset of chains: ",
                        paste(head(bad_variables, 20), collapse = ", "),
                        if (length(bad_variables) > 20) " ..." else ""
                    ),
                    call. = FALSE
                )
            }

            if (any(globally_constant)) {
                message(
                    sum(globally_constant),
                    " globally constant variables were excluded from ",
                    "the autocorrelation check."
                )
            }

            variables_to_check <- which(!globally_constant)

            mean_rho1 <- rep(NA_real_, n_variables)
            names(mean_rho1) <- variable_names

            if (length(variables_to_check) == 0L) {
                return(mean_rho1)
            }

            rho1_by_chain <- matrix(
                NA_real_,
                nrow = n_chains,
                ncol = length(variables_to_check),
                dimnames = list(
                    chain = seq_len(n_chains),
                    variable = variable_names[variables_to_check]
                )
            )

            for (j in seq_along(variables_to_check)) {
                variable_index <- variables_to_check[j]

                for (chain_index in seq_len(n_chains)) {
                    z <- x[, chain_index, variable_index]

                    rho1_by_chain[chain_index, j] <- stats::acf(
                        z,
                        lag.max = 1,
                        plot = FALSE,
                        demean = TRUE
                    )$acf[2]
                }
            }

            mean_rho1[variables_to_check] <- colMeans(abs(rho1_by_chain))

            checked_rho1 <- mean_rho1[variables_to_check]
            if (anyNA(checked_rho1) || any(!is.finite(checked_rho1))) {
                stop(
                    "Autocorrelation was unexpectedly undefined for a nonconstant variable.",
                    call. = FALSE
                )
            }

            mean_rho1
        },
        check_rhat = function(draws, rhat, threshold = 1.01) {
            draws <- posterior::as_draws_array(draws)

            variable_names <- dimnames(draws)$variable

            if (length(rhat) != length(variable_names)) {
                stop(
                    "R-hat values do not match the number of draw variables.",
                    call. = FALSE
                )
            }
            if (!is.null(names(rhat))) {
                missing_rhat <- setdiff(variable_names, names(rhat))
                if (length(missing_rhat) > 0L) {
                    stop(
                        paste0(
                            "R-hat values are missing for: ",
                            paste(missing_rhat, collapse = ", ")
                        ),
                        call. = FALSE
                    )
                }
                rhat <- rhat[variable_names]
            }

            all_finite <- apply(
                draws,
                3,
                function(z) all(is.finite(z))
            )

            if (any(!all_finite)) {
                stop(
                    "Non-finite posterior draws were found.",
                    call. = FALSE
                )
            }

            globally_constant <- apply(
                draws,
                3,
                self$is_constant
            )

            constant_by_chain <- apply(
                draws,
                c(2, 3),
                self$is_constant
            )

            partially_constant <- apply(
                constant_by_chain,
                2,
                any
            ) &
                !globally_constant

            if (any(partially_constant)) {
                stop(
                    paste0(
                        "Variables were constant in only some chains, or ",
                        "constant at different values across chains: ",
                        paste(
                            variable_names[partially_constant],
                            collapse = ", "
                        )
                    ),
                    call. = FALSE
                )
            }

            # names(rhat) <- variable_names

            variables_to_check <- !globally_constant

            checked_rhat <- rhat[variables_to_check]

            if (
                anyNA(checked_rhat) ||
                    any(!is.finite(checked_rhat))
            ) {
                stop(
                    "R-hat was unexpectedly undefined for a nonconstant variable.",
                    call. = FALSE
                )
            }

            all(checked_rhat < threshold)
        },
        get_diagnostics = function(stan_fit, to_keep = NULL) {
            diag_summ <- rstan::get_sampler_params(stan_fit, inc_warmup = FALSE)
            draws <- posterior::as_draws_df(stan_fit)
            if (is.null(to_keep)) {
                to_keep <- self$get_posterior_dims(stan_fit@model_name)
            }
            if (length(to_keep) != 0) {
                draws <- posterior::subset_draws(draws, names(to_keep))
            } else {
                warning(paste0(
                    "Including all variables in diagnostics. ",
                    "This likely means that the posterior object ",
                    "does not contain necessary model dimensions."
                ))
            }
            summ <- posterior::summarize_draws(draws)
            diagnostic_names <- summ$variable
            diagnostics <- list(
                diagnostic_information = list(names = diagnostic_names),
                ndraws = posterior::ndraws(draws),
                nchains = posterior::nchains(draws),
                effective_sample_size_bulk = stats::setNames(
                    summ$ess_bulk,
                    diagnostic_names
                ),
                effective_sample_size_tail = stats::setNames(
                    summ$ess_tail,
                    diagnostic_names
                ),
                rhat = stats::setNames(summ$rhat, diagnostic_names),
                divergent_transitions = sapply(
                    diag_summ,
                    function(x) {
                        sum(x[, "divergent__"])
                    }
                ),
                efmi = rstan::get_bfmi(stan_fit),
                mean_lag1_ac = self$compute_mean_lag1_ac(draws),
                mean_lag1_ac_posterior = self$compute_ac(draws)
            )
            diagnostics
        },
        get_posterior_json = function(
            posterior_name,
            local_posterior_path = "/posterior_database/posteriors/"
        ) {
            posterior_name <- paste0(posterior_name, ".json")
            posterior_path <- normalizePath(file.path(
                self$get_pdb_path(),
                sub("^[/\\\\]+", "", local_posterior_path),
                posterior_name
            ))
            jsonlite::read_json(posterior_path, simplifyVector = TRUE)
        },
        get_posterior_dims = function(po) {
            self$get_posterior_json(po)$dimensions
        },
        is_constant = function(z, tolerance = sqrt(.Machine$double.eps)) {
            if (!all(is.finite(z))) {
                return(FALSE)
            }

            scale <- max(1, max(abs(z)))
            diff(range(z)) <= tolerance * scale
        },
        get_reference_info = function(x) {
            attr(x, "info", exact = TRUE)
        },
        set_reference_info = function(x, value) {
            if (!is.list(value)) {
                stop("Reference information must be a list.", call. = FALSE)
            }
            attr(x, "info") <- value
            x
        },
        search_data = function(query) {
            private$search_database_entries(
                query,
                "posterior_database/data/data",
                "\\.json(?:\\.zip)?$",
                "data"
            )
        },
        search_data_keywords = function(query) {
            private$search_info_keywords(
                query,
                "posterior_database/data/info",
                "data"
            )
        },
        search_model = function(query) {
            private$search_database_entries(
                query,
                "posterior_database/models/info",
                "\\.info\\.json$",
                "model"
            )
        },
        search_model_keywords = function(query) {
            private$search_info_keywords(
                query,
                "posterior_database/models/info",
                "model"
            )
        },
        search_posterior = function(query) {
            private$search_database_entries(
                query,
                "posterior_database/posteriors",
                "\\.json$",
                "posterior"
            )
        },
        search_reference_draws = function(query) {
            private$search_database_entries(
                query,
                paste0(
                    "posterior_database/reference_posteriors/",
                    "draws/draws"
                ),
                "\\.json\\.zip$",
                "reference-draw"
            )
        },
        add_bibtex_entry = function(bibtex_str) {
            checkmate::assert_string(bibtex_str, min.chars = 1L)

            extract_keys <- function(x) {
                pattern <- paste0(
                    "@([[:alpha:]]+)\\s*[\\{(]\\s*",
                    "([^,[:space:]]+)\\s*,"
                )
                matches <- regmatches(
                    x,
                    gregexpr(pattern, x, perl = TRUE, ignore.case = TRUE)
                )[[1L]]
                if (identical(matches, character(0))) {
                    return(character())
                }
                types <- sub(
                    pattern,
                    "\\1",
                    matches,
                    perl = TRUE,
                    ignore.case = TRUE
                )
                keys <- sub(
                    pattern,
                    "\\2",
                    matches,
                    perl = TRUE,
                    ignore.case = TRUE
                )
                keys[!tolower(types) %in% c("comment", "preamble", "string")]
            }

            new_keys <- extract_keys(bibtex_str)
            if (length(new_keys) != 1L) {
                stop(
                    "`bibtex_str` must contain exactly one keyed BibTeX entry.",
                    call. = FALSE
                )
            }
            reference_path <- private$get_reference_path()
            existing_text <- if (file.exists(reference_path)) {
                paste(readLines(reference_path, warn = FALSE), collapse = "\n")
            } else {
                ""
            }
            if (tolower(new_keys) %in% tolower(extract_keys(existing_text))) {
                message(
                    "BibTeX entry `",
                    new_keys,
                    "` already exists; nothing was written."
                )
                return(invisible(FALSE))
            }

            write(
                paste0("\n\n", bibtex_str, "\n"),
                file = reference_path,
                append = TRUE
            )
            self$refresh()
            invisible(TRUE)
        },
        set_data = function(d) {
            self$data <- d
        },
        get_data = function() {
            self$data
        },
        set_model_code = function(c) {
            self$model_code <- c
        },
        get_model_code = function() {
            self$model_code
        },
        get_stan_model_code_path = function(
            local_model_path = "/posterior_database/models/stan/"
        ) {
            normalizePath(file.path(
                self$get_pdb_path(),
                sub("^[/\\\\]+", "", local_model_path)
            ))
        },
        get_data_path = function(
            local_data_path = "/posterior_database/data/data/"
        ) {
            normalizePath(file.path(
                self$get_pdb_path(),
                sub("^[/\\\\]+", "", local_data_path)
            ))
        },
        get_posterior_modeldata = function(posterior_name) {
            md_files <- private$get_posterior_modeldata_files(posterior_name)
            files_exist <- vapply(md_files, file.exists, logical(1))
            if (!all(files_exist)) {
                stop(
                    paste0(
                        "The data or model file does not exist. ",
                        "Please ensure you specified the correct posterior name."
                    ),
                    call. = FALSE
                )
            }
            posterior_model <- rstan::stan_model(
                file = md_files$model_file,
                model_name = posterior_name
            )
            posterior_model_code <- paste(
                readLines(md_files$model_file),
                collapse = "\n"
            )
            posterior_data <- private$copy_to_tempdir(md_files$data_file)
            list(
                data = posterior_data,
                stan_model = posterior_model,
                model_code = posterior_model_code,
                stan_file = md_files$model_file
            )
        },
        # get_reference_posterior_info = function(
        #     po,
        #     local_posterior_info_path = "/posterior_database/reference_posteriors/draws/info/"
        # ) {
        #     rpi_path <- normalizePath(file.path(self$get_pdb_path(), )
        # },
        compute_reference_summary_statistics = function(stan_fit) {
            stan_info <- self$get_reference_info(stan_fit)
            if (is.null(stan_info) || is.null(stan_info$checks_made)) {
                stop(
                    "Draws must be checked before computing summaries.",
                    call. = FALSE
                )
            }
            required_checks <- c(
                "ndraws_is_10k",
                "nchains_is_gte_4",
                "r_hat_below_1_01",
                "efmi_above_0_2",
                "abs_mean_lag1_ac_below_0_05"
            )
            failed_checks <- required_checks[
                !vapply(
                    stan_info$checks_made[required_checks],
                    isTRUE,
                    logical(1)
                )
            ]
            if (length(failed_checks) > 0L) {
                stop(
                    "Required checks failed: ",
                    paste(failed_checks, collapse = ", "),
                    call. = FALSE
                )
            }
            draws <- posterior::subset_draws(
                posterior::as_draws_array(stan_fit),
                variable = names(self$get_posterior_dims(stan_fit@model_name))
            )
            mean_summary <- posterior::summarize_draws(
                draws,
                "mean",
                "mcse_mean"
            )
            squared_draws <- draws
            squared_draws[] <- squared_draws[]^2
            squared_summary <- posterior::summarize_draws(
                squared_draws,
                "mean",
                "mcse_mean"
            )
            list(
                mean_value = list(
                    names = as.character(mean_summary$variable),
                    mean_value = as.numeric(mean_summary$mean),
                    mcse_mean = as.numeric(mean_summary$mcse_mean)
                ),
                mean_squared_value = list(
                    names = as.character(squared_summary$variable),
                    mean_squared_value = as.numeric(squared_summary$mean),
                    mcse_mean = as.numeric(squared_summary$mcse_mean)
                )
            )
        },
        write_summary_statistics_from_stan_fit = function(
            stan_fit,
            overwrite = FALSE,
            verify = TRUE
        ) {
            checkmate::assert_flag(overwrite)
            checkmate::assert_flag(verify)
            summaries <- self$compute_reference_summary_statistics(stan_fit)
            summary_info <- self$get_reference_info(stan_fit)
            summary_info$versions$r_summary_statistic <- paste0(
                "posterior R package, version ",
                utils::packageVersion("posterior")
            )
            paths <- lapply(names(summaries), function(type) {
                private$get_summary_statistic_paths(stan_fit@model_name, type)
            })
            names(paths) <- names(summaries)
            destinations <- unlist(paths, use.names = FALSE)
            if (!overwrite && any(file.exists(destinations))) {
                stop(
                    "One or more summary-statistic files already exist.",
                    call. = FALSE
                )
            }
            write_one <- function(x, path) {
                dir.create(
                    dirname(path),
                    recursive = TRUE,
                    showWarnings = FALSE
                )
                jsonlite::write_json(
                    x,
                    path,
                    pretty = TRUE,
                    auto_unbox = TRUE,
                    null = "null",
                    digits = NA
                )
            }
            for (type in names(summaries)) {
                write_one(summary_info, paths[[type]]$info)
                write_one(summaries[[type]], paths[[type]]$value)
                if (verify) {
                    value <- jsonlite::read_json(
                        paths[[type]]$value,
                        simplifyVector = TRUE
                    )
                    info <- jsonlite::read_json(
                        paths[[type]]$info,
                        simplifyVector = TRUE
                    )
                    if (
                        !identical(info$name, summary_info$name) ||
                            !identical(value$names, summaries[[type]]$names) ||
                            !isTRUE(all.equal(
                                value[[type]],
                                summaries[[type]][[type]],
                                check.attributes = FALSE
                            )) ||
                            !isTRUE(all.equal(
                                value$mcse_mean,
                                summaries[[type]]$mcse_mean,
                                check.attributes = FALSE
                            ))
                    ) {
                        stop(
                            "Summary verification failed for ",
                            type,
                            ".",
                            call. = FALSE
                        )
                    }
                }
            }
            invisible(paths)
        },
        get_rpi_path = function(
            rpi_name,
            local_rpi_path = "/posterior_database/reference_posteriors/draws/info"
        ) {
            rpi_name <- paste0(rpi_name, ".info.json")
            normalizePath(
                file.path(
                    self$get_pdb_path(),
                    sub("^[/\\\\]+", "", local_rpi_path),
                    rpi_name
                ),
                mustWork = FALSE
            )
        },
        get_rp_path = function(
            rp_name,
            local_rp_path = "/posterior_database/reference_posteriors/draws/draws"
        ) {
            rp_name <- paste0(rp_name, ".json.zip")
            normalizePath(
                file.path(
                    self$get_pdb_path(),
                    sub("^[/\\\\]+", "", local_rp_path),
                    rp_name
                ),
                mustWork = FALSE
            )
        },
        set_stan_file = function(sf) {
            self$stan_file <- sf
        },
        get_stan_file = function() {
            self$stan_file
        },
        set_posterior = function(p) {
            self$posterior <- p
        },
        get_posterior = function() {
            self$posterior
        },
        set_rp = function(rp) {
            self$rp <- rp
        },
        get_rp = function() {
            self$rp
        },
        get_added_by = function() {
            self$adder
        },
        set_added_by = function(added_by) {
            self$adder <- added_by
        },
        refresh = function() {
            private$pdb <- posteriordb::pdb_local(path = self$path)
            invisible(self)
        },
        get_pdb = function() {
            private$pdb
        },
        get_pdb_path = function() {
            self$path
        }
    ),

    private = list(
        pdb = NULL,
        get_posterior_modeldata_files = function(posterior_name) {
            data_model_name <- strsplit(posterior_name, split = "-")[[1]]
            if (length(data_model_name) < 2L) {
                stop(
                    "`posterior_name` must identify both data and model.",
                    call. = FALSE
                )
            }
            data_file_name <- paste0(data_model_name[1], ".json.zip")
            model_file_name <- paste0(data_model_name[2], ".stan")
            list(
                "model_file" = normalizePath(file.path(
                    self$get_stan_model_code_path(),
                    model_file_name
                )),
                "data_file" = normalizePath(file.path(
                    self$get_data_path(),
                    data_file_name
                ))
            )
        },
        get_summary_statistic_paths = function(posterior_name, type) {
            supported <- c("mean_value", "mean_squared_value")
            checkmate::assert_choice(type, supported)
            root <- file.path(
                self$get_pdb_path(),
                "posterior_database",
                "reference_posteriors",
                "summary_statistics",
                type
            )
            list(
                info = file.path(
                    root,
                    "info",
                    paste0(posterior_name, ".info.json")
                ),
                value = file.path(
                    root,
                    type,
                    paste0(posterior_name, ".json")
                )
            )
        },
        copy_to_tempdir = function(
            file_path,
            return_obj = TRUE,
            overwrite = TRUE
        ) {
            is_zip <- grepl("\\.zip$", file_path, ignore.case = TRUE)
            task_dir <- tempfile(pattern = "pdb-entry-")
            dir.create(task_dir)
            on.exit(
                unlink(task_dir, recursive = TRUE, force = TRUE),
                add = TRUE
            )
            copied_path <- file.path(task_dir, basename(file_path))
            copied <- file.copy(
                from = file_path,
                to = copied_path,
                overwrite = overwrite
            )
            if (!copied) {
                stop(
                    "Failed to copy file to a temporary directory.",
                    call. = FALSE
                )
            }
            if (return_obj && is_zip) {
                archive <- utils::unzip(copied_path, list = TRUE)
                if (nrow(archive) != 1L) {
                    stop(
                        "Expected a ZIP archive containing one file.",
                        call. = FALSE
                    )
                }
                extracted <- utils::unzip(
                    copied_path,
                    files = archive$Name[[1]],
                    exdir = task_dir
                )
                jsonlite::read_json(extracted, simplifyVector = TRUE)
            } else if (return_obj) {
                jsonlite::read_json(copied_path, simplifyVector = TRUE)
            } else {
                destination <- tempfile(
                    pattern = "pdb-entry-copy-",
                    fileext = paste0(".", tools::file_ext(file_path))
                )
                if (!file.copy(copied_path, destination, overwrite = TRUE)) {
                    stop("Failed to retain the temporary copy.", call. = FALSE)
                }
                destination
            }
        },
        search_database_entries = function(
            query,
            directory,
            suffix_pattern,
            entry_type
        ) {
            checkmate::assert_string(query, min.chars = 1L)
            search_path <- file.path(self$path, directory)
            if (!dir.exists(search_path)) {
                stop(
                    "Could not find the ",
                    entry_type,
                    " directory: ",
                    search_path,
                    call. = FALSE
                )
            }

            files <- list.files(search_path, full.names = FALSE)
            files <- files[grepl(
                suffix_pattern,
                files,
                ignore.case = TRUE,
                perl = TRUE
            )]
            names <- sub(
                suffix_pattern,
                "",
                files,
                ignore.case = TRUE,
                perl = TRUE
            )
            matches <- sort(unique(names[
                grepl(tolower(query), tolower(names), fixed = TRUE)
            ]))

            if (length(matches) == 0L) {
                message("No ", entry_type, " names matched `", query, "`.")
            }
            matches
        },
        search_info_keywords = function(query, directory, entry_type) {
            checkmate::assert_string(query, min.chars = 1L)
            search_path <- file.path(self$path, directory)
            if (!dir.exists(search_path)) {
                stop(
                    "Could not find the ",
                    entry_type,
                    " info directory: ",
                    search_path,
                    call. = FALSE
                )
            }

            info_files <- list.files(
                search_path,
                pattern = "\\.info\\.json$",
                full.names = TRUE,
                ignore.case = TRUE
            )
            matches <- vapply(
                info_files,
                function(info_file) {
                    info <- tryCatch(
                        jsonlite::fromJSON(info_file, simplifyVector = FALSE),
                        error = function(e) NULL
                    )
                    if (is.null(info)) {
                        return(NA_character_)
                    }
                    keywords <- as.character(unlist(
                        info$keywords,
                        recursive = TRUE,
                        use.names = FALSE
                    ))
                    if (
                        length(keywords) == 0L ||
                            !any(grepl(
                                tolower(query),
                                tolower(keywords),
                                fixed = TRUE
                            ))
                    ) {
                        return(NA_character_)
                    }
                    as.character(info$name)[1L]
                },
                character(1L),
                USE.NAMES = FALSE
            )
            matches <- sort(unique(matches[!is.na(matches)]))

            if (length(matches) == 0L) {
                message(
                    "No ",
                    entry_type,
                    " keywords matched `",
                    query,
                    "`."
                )
            }
            matches
        },
        get_sampling_version_info = function() {
            M <- file.path(
                Sys.getenv("HOME"),
                ".R",
                ifelse(
                    .Platform$OS.type == "windows",
                    "Makevars.win",
                    "Makevars"
                )
            )
            Mfile <- if (file.exists(M)) {
                paste(readLines(M), collapse = "\n")
            } else {
                "[Could not find Makevar file]"
            }
            list(
                rstan_version = paste("rstan", utils::packageVersion("rstan")),
                r_Makevars = paste(Mfile, collapse = "\n"),
                r_version = R.version$version.string,
                r_session = paste(
                    utils::capture.output(print(utils::sessionInfo())),
                    collapse = "\n"
                )
            )
        },
        get_reference_path = function(
            local_bib_path = "posterior_database/bibliography/references.bib"
        ) {
            file.path(self$path, local_bib_path)
        }
    )
)
