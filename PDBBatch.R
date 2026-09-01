library(R6)

PDBBatch <- R6Class(
    "PDBBatch",
    public = list(
        path = NULL,
        added_by = NULL,
        initialize = function(path, added_by = "Gerald Press") {
            req <- c(
                "PDBEntryBuilder",
                "PDBDataEntry",
                "PDBModelEntry",
                "PDBPosteriorEntry"
            )
            miss <- req[!vapply(req, exists, logical(1), inherits = TRUE)]
            if (length(miss)) {
                stop(
                    "Source PDBEntryBuilder.R and PDBEntry.R first. Missing: ",
                    paste(miss, collapse = ", "),
                    call. = FALSE
                )
            }
            self$path <- normalizePath(path, mustWork = TRUE)
            self$added_by <- added_by
            private$builder <- PDBEntryBuilder$new(self$path, added_by)
            invisible(self)
        },
        add_data = function(
            name,
            data = NULL,
            info = NULL,
            raw_data = data,
            data_info = info
        ) {
            private$unique(name, private$data, "data")
            x <- PDBDataEntry$new(name, raw_data, data_info)
            private$data[[name]] <- x
            invisible(x)
        },
        add_model = function(name, stan_file, info = NULL, model_info = info) {
            private$unique(name, private$models, "model")
            x <- PDBModelEntry$new(name, stan_file, model_info)
            private$models[[name]] <- x
            invisible(x)
        },
        add_posterior = function(
            data_name,
            model_name,
            name = paste(data_name, model_name, sep = "-"),
            keywords = NULL,
            references = NULL,
            framework = "stan",
            dimensions = NULL,
            added_date = Sys.Date(),
            added_by = self$added_by,
            parameters_include = NULL,
            parameters_exclude = NULL,
            reference_sampling = NULL
        ) {
            private$known(data_name, private$data, "data")
            private$known(model_name, private$models, "model")
            private$unique(name, private$posteriors, "posterior")
            spec <- list(
                keywords = keywords,
                references = references,
                framework = framework,
                dimensions = dimensions,
                added_date = added_date,
                added_by = added_by,
                parameters_include = parameters_include,
                parameters_exclude = parameters_exclude
            )
            x <- PDBPosteriorEntry$new(
                name,
                private$data[[data_name]],
                private$models[[model_name]],
                spec,
                reference_sampling
            )
            private$posteriors[[name]] <- x
            invisible(x)
        },
        add_entry = function(...) self$add_posterior(...),
        add_grid = function(
            data_names = NULL,
            model_names = NULL,
            common = list(),
            overrides = list()
        ) {
            if (is.null(data_names)) {
                data_names <- names(private$data)
            }
            if (is.null(model_names)) {
                model_names <- names(private$models)
            }
            if (!is.list(common) || !is.list(overrides)) {
                stop("`common` and `overrides` must be lists.", call. = FALSE)
            }
            private$known(data_names, private$data, "data")
            private$known(model_names, private$models, "model")
            if (!length(data_names) || !length(model_names)) {
                stop("Grid requires data and models.", call. = FALSE)
            }
            g <- expand.grid(
                data_name = data_names,
                model_name = model_names,
                stringsAsFactors = FALSE
            )
            n <- paste(g$data_name, g$model_name, sep = "-")
            u <- setdiff(names(overrides), n)
            if (length(u)) {
                stop(
                    "Unknown overrides: ",
                    paste(u, collapse = ", "),
                    call. = FALSE
                )
            }
            d <- intersect(n, names(private$posteriors))
            if (length(d)) {
                stop(
                    "Existing posteriors: ",
                    paste(d, collapse = ", "),
                    call. = FALSE
                )
            }
            allowed <- names(formals(self$add_posterior))
            private$settings(common, allowed, "common")
            for (k in names(overrides)) {
                if (!is.list(overrides[[k]])) {
                    stop("Overrides must be lists.", call. = FALSE)
                }
                private$settings(overrides[[k]], allowed, k)
            }
            out <- setNames(vector("list", nrow(g)), n)
            for (i in seq_len(nrow(g))) {
                o <- overrides[[n[[i]]]]
                if (is.null(o)) {
                    o <- list()
                }
                a <- utils::modifyList(common, o)
                a$data_name <- g$data_name[[i]]
                a$model_name <- g$model_name[[i]]
                a$name <- n[[i]]
                out[[i]] <- do.call(self$add_posterior, a)
            }
            invisible(out)
        },
        validate = function(deep = FALSE) {
            private$flag(deep, "deep")
            for (x in private$data) {
                private$validate_data(x, deep)
            }
            for (x in private$models) {
                private$validate_model(x, deep)
            }
            for (x in private$posteriors) {
                private$validate_posterior(x, deep)
            }
            invisible(self$results())
        },
        write_all = function(
            overwrite = FALSE,
            on_error = c("continue", "stop")
        ) {
            private$flag(overwrite, "overwrite")
            on_error <- match.arg(on_error)
            self$validate(FALSE)
            for (x in private$data) {
                if (x$results$validation$status != "valid") {
                    next
                }
                ok <- private$attempt(x, "writing", function() {
                    private$builder$add_data(x$raw_data, x$data_info, overwrite)
                })
                if (!ok && on_error == "stop") return(invisible(self$results()))
            }
            for (x in private$models) {
                if (x$results$validation$status != "valid") {
                    next
                }
                ok <- private$attempt(x, "writing", function() {
                    private$builder$add_model_code(
                        x$stan_file,
                        x$model_info,
                        overwrite
                    )
                })
                if (!ok && on_error == "stop") return(invisible(self$results()))
            }
            for (x in private$posteriors) {
                if (
                    x$results$validation$status != "valid" ||
                        is.null(x$data_entry$pdb_data) ||
                        is.null(x$model_entry$pdb_model_code)
                ) {
                    x$set_result(
                        "writing",
                        "skipped",
                        error = "A dependency was not written."
                    )
                    next
                }
                ok <- private$attempt(x, "writing", function() {
                    private$builder$add_posterior(
                        private$spec(x),
                        overwrite,
                        FALSE
                    )
                })
                if (!ok && on_error == "stop") return(invisible(self$results()))
            }
            invisible(self$results())
        },
        compute_references = function(
            entries = NULL,
            sampling_args = NULL,
            check_and_write = TRUE,
            overwrite = FALSE,
            on_error = c("continue", "stop"),
            ...
        ) {
            on_error <- match.arg(on_error)
            if (is.null(entries)) {
                entries <- names(private$posteriors)
            }
            private$known(entries, private$posteriors, "posterior")
            extra <- list(...)
            for (n in entries) {
                x <- private$posteriors[[n]]
                inf <- sampling_args
                if (is.null(inf)) {
                    inf <- x$reference_sampling
                }
                ri <- private$reference_info(x, inf, extra)
                x$set_reference_info(ri)
                ok <- private$attempt(x, "reference", function() {
                    do.call(
                        private$builder$compute_reference,
                        c(
                            list(
                                posterior_name = x$name,
                                sampling_args = inf,
                                check_and_write = check_and_write,
                                overwrite = overwrite
                            ),
                            extra
                        )
                    )
                })
                if (!ok && on_error == "stop") break
            }
            invisible(self$results())
        },
        results = function() {
            f <- function(xs, type) {
                unlist(
                    lapply(xs, function(x) {
                        lapply(names(x$results), function(s) {
                            c(
                                list(type = type, name = x$name, stage = s),
                                x$results[[s]]
                            )
                        })
                    }),
                    recursive = FALSE
                )
            }
            c(
                f(private$data, "data"),
                f(private$models, "model"),
                f(private$posteriors, "posterior")
            )
        },
        get_data = function(name = NULL) {
            private$get(private$data, name, "data")
        },
        get_model = function(name = NULL) {
            private$get(private$models, name, "model")
        },
        get_posterior = function(name = NULL) {
            private$get(private$posteriors, name, "posterior")
        },
        get_entry = function(name = NULL) self$get_posterior(name),
        get_builder = function() private$builder
    ),
    private = list(
        builder = NULL,
        data = list(),
        models = list(),
        posteriors = list(),
        unique = function(n, r, t) {
            .pdb_name(n)
            if (n %in% names(r)) stop("Duplicate ", t, ": ", n, call. = FALSE)
        },
        known = function(n, r, t) {
            u <- setdiff(n, names(r))
            if (length(u)) {
                stop(
                    "Unknown ",
                    t,
                    ": ",
                    paste(u, collapse = ", "),
                    call. = FALSE
                )
            }
        },
        flag = function(x, n) {
            if (!is.logical(x) || length(x) != 1L || is.na(x)) {
                stop("`", n, "` must be TRUE or FALSE.", call. = FALSE)
            }
        },
        settings = function(x, a, n) {
            u <- setdiff(names(x), a)
            if (length(u)) {
                stop(
                    "Unknown settings for ",
                    n,
                    ": ",
                    paste(u, collapse = ", "),
                    call. = FALSE
                )
            }
        },
        get = function(r, n, t) {
            if (is.null(n)) {
                return(r)
            }
            private$known(n, r, t)
            r[[n]]
        },
        attempt = function(x, s, f) {
            tryCatch(
                {
                    o <- f()
                    x$set_result(s, "added", o)
                    TRUE
                },
                error = function(e) {
                    x$set_result(s, "failed", error = e)
                    FALSE
                }
            )
        },
        validate_data = function(x, deep) {
            e <- tryCatch(
                {
                    n <- names(x$raw_data)
                    if (
                        is.null(n) ||
                            anyNA(n) ||
                            any(n == "") ||
                            anyDuplicated(n)
                    ) {
                        stop("Raw data must have unique non-empty names.")
                    }
                    if (deep) {
                        o <- private$builder$create_data(
                            x$raw_data,
                            x$data_info
                        )
                        x$set_result("creation", "added", o)
                    }
                    NULL
                },
                error = identity
            )
            x$set_result(
                "validation",
                if (is.null(e)) "valid" else "invalid",
                error = e
            )
        },
        validate_model = function(x, deep) {
            e <- tryCatch(
                {
                    if (!file.exists(x$stan_file)) {
                        stop("Stan file does not exist: ", x$stan_file)
                    }
                    if (deep) {
                        o <- private$builder$create_model_code(
                            x$stan_file,
                            x$model_info
                        )
                        x$set_result("creation", "added", o)
                    }
                    NULL
                },
                error = identity
            )
            x$set_result(
                "validation",
                if (is.null(e)) "valid" else "invalid",
                error = e
            )
        },
        validate_posterior = function(x, deep) {
            e <- tryCatch(
                {
                    i <- x$posterior_spec$parameters_include
                    e <- x$posterior_spec$parameters_exclude
                    c <- intersect(i, e)
                    if (length(c)) {
                        stop(
                            "Included and excluded: ",
                            paste(c, collapse = ", ")
                        )
                    }
                    if ("lp__" %in% i) {
                        stop("`lp__` cannot be included.")
                    }
                    if (deep) {
                        if (
                            is.null(x$data_entry$pdb_data) ||
                                is.null(x$model_entry$pdb_model_code)
                        ) {
                            stop("Invalid dependency.")
                        }
                        o <- private$builder$add_posterior(
                            private$spec(x),
                            dry_run = TRUE
                        )
                        private$builder$check_posterior(o, TRUE)
                        x$set_result("creation", "added", o)
                    }
                    NULL
                },
                error = identity
            )
            x$set_result(
                "validation",
                if (is.null(e)) "valid" else "invalid",
                error = e
            )
        },
        spec = function(x) {
            c(
                list(
                    pdb_data = x$data_entry$pdb_data,
                    pdb_model_code = x$model_entry$pdb_model_code
                ),
                x$posterior_spec
            )
        },
        reference_info = function(x, inf, extra) {
            if (is.null(inf)) {
                inf <- list(
                    method = "stan_sampling",
                    method_arguments = list(
                        chains = 10,
                        iter = 20000,
                        warmup = 10000,
                        thin = 10,
                        refresh = 10000,
                        seed = if (is.null(extra$seed)) 123 else extra$seed,
                        control = if (is.null(extra$control_args)) {
                            list(adapt_delta = .9)
                        } else {
                            extra$control_args
                        }
                    )
                )
            }
            z <- list(
                name = x$name,
                inference = inf,
                diagnostics = NULL,
                checks_made = NULL,
                comments = extra$comments,
                added_by = self$added_by,
                added_date = Sys.Date(),
                versions = NULL
            )
            if (requireNamespace("posteriordb", quietly = TRUE)) {
                tryCatch(
                    posteriordb::as.pdb_reference_posterior_info(z),
                    error = function(e) z
                )
            } else {
                z
            }
        }
    )
)
