library(R6)

PDBEntryBuilder <- R6Class(
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
            require(rstan)
            require(posterior)
            rstan_options(auto_write = auto_write)
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
        #TO DO
        compute_reference = function(
            posterior_name,
            sampling_args = NULL,
            seed = 123,
            control_args = list(adapt_delta = 0.9),
            comments = NULL,
            added_by = self$adder,
            check_and_write = TRUE,
            overwrite = FALSE
        ) {
            if (is.null(sampling_args)) {
                ri <- list(
                    name = posterior_name,
                    inference = list(
                        method = "stan_sampling",
                        method_arguments = list(
                            chains = 10,
                            iter = 20000,
                            warmup = 10000,
                            thin = 10,
                            refresh = 10000,
                            seed = seed,
                            control = control_args
                        ),
                        diagnostics = NULL,
                        checks_made = NULL,
                        comments = comments,
                        added_by = self$adder,
                        added_date = Sys.Date(),
                        versions = NULL
                    )
                )
            } else {
                ri <- list(
                    name = posterior_name,
                    inference = sampling_args,
                    diagnostics = NULL,
                    checks_made = NULL,
                    comments = comments,
                    added_by = self$adder,
                    added_date = Sys.Date(),
                    versions = NULL
                )
            }
            rpi <- posteriordb::as.pdb_reference_posterior_info(ri)
            rp <- posteriordb::compute_reference_posterior_draws(
                rpi,
                private$pdb
            )
            if (check_and_write) {
                rp <- self$check_reference_draws(rp)
                posteriordb::write_pdb(rp, private$pdb, overwrite = overwrite)
            }
            # Uses private$pdb
            invisible(rp)
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
            info(stan_fit) <- list(
                name = posterior_name,
                inference = list(
                    method = "stan_sampling",
                    method_arguments = inference
                ),
                diagnostics = self$get_diagnostics(stan_fit, to_keep = ),
                checks_made = list(),
                comments = comments,
                added_by = added_by,
                added_date = Sys.Date(),
                versions = private$get_sampling_version_info()
            )
            if (auto_check) {
                stan_fit <- self$check_draws_from_stanfit(stan_fit)
            }
            #TO DO: Write_reference_draws_from_stan_fit
            # if (write) {
            #     self$write_reference_draws(stan_fit, overwrite = overwrite)
            # }
            stan_fit
        },
        write_rpi_from_stan_fit = function(stan_fit) {
            if (is.null(info(stan_fit))) {
                stop(
                    "`stan_fit` object must contain reference draw info. Call `add_rpi_from_stanfit()`",
                    call. = FALSE
                )
            }
            info <- info(stan_fit)
            if (is.null(info$checks_made)) {
                stop(
                    "Draws must be checked before writing. Call `check_draws_from_stanfit()`",
                    call. = FALSE
                )
            }
            checks <- unlist(info$checks_made)
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
            jsonlite::write_json(
                info(stan_fit),
                self$get_rpi_path(stan_fit@model_name)
            )
        },
        write_rpd_from_stan_fit = function(stan_fit) {
            if (!file.exists(self$get_rpi_path(stan_fit@model_name))) {
                stop(
                    "Make sure to write the reference posterior information to disk before writing the draws",
                    call. = FALSE
                )
            }
            to_keep <- self$get_posterior_dims(stan_fit@model_name)
            draws <- posterior::subset_draws(
                posterior::as_draws_array(stan_fit),
                to_keep
            )
            rp_path <- self$get_rp_path(stan_fit@model_name)
            jsonlite::write_json(draws, rp_path)
            zip(rp_path)
            file.remove(rp_path)
            invisible(TRUE)
        },
        add_reference_draw_info_stanfit = function(stan_fit) {},
        check_reference_draws = function(rp) {
            if (is.null(rp)) {
                if (is.null(self$get_rp())) {
                    stop("Reference draws not found.", call. = FALSE)
                }
                posteriordb::check_reference_posterior_draws(x = self$get_rp())
            } else {
                posteriordb::check_reference_posterior_draws(x = rp)
            }
        },
        check_draws_from_stanfit = function(stan_fit, attach = FALSE) {
            if (is.null(info(stan_fit)$diagnostics)) {
                diagnostics <- self$get_diagnostics(stan_fit)
                info(stan_fit)$diagnostics <- self$get_diagnostics(stan_fit)
            } else {
                diagnostics <- info(stan_fit)$diagnostics
            }
            if (any(diagnostics$divergent_transitions > 0)) {
                stop(
                    "There were divergent transitions during sampling.",
                    call. = FALSE
                )
            }
            ess_failures <- self$get_ess_bounds_failures(diagnostics)
            ess_within_bounds <- ess_failures$total_count == 0L
            checks_made <- list(
                ndraws_is_10k = diagnostics$ndraws == 10000,
                nchains_is_gte_4 = diagnostics$nchains >= 4,
                r_hat_below_1_01 = self$check_rhat(
                    stan_fit,
                    diagnostics$rhat
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
            if (attach) {
                info(stan_fit)$checks_made <- checks_made
            } else {
                invisible(checks_made)
            }
        },
        check_draws = function(rp, write_mean_ac = FALSE) {
            rpi <- info(rp)
            if (is.null(rpi)) {
                stop(
                    "Posterior Reference Draws do not have required information attached.",
                    call. = FALSE
                )
            }
            if (inherits(rpi, "list")) {
                message(
                    "Attempting to convert information list into `pdb_reference_posterior_info` object..."
                )
                rpi <- as.pdb_reference_posterior_info(rpi)
            }
            if (!inherits(rpi, "pdb_reference_posterior_info")) {
                stop(
                    "Posterior information needs to be embedded as a `pdb_reference_posterior_info` object.",
                    call. = FALSE
                )
            }
            checks_made <- rpi$checks_made
            if (is.null(rpi$diagnostics)) {
                rpi$diagnostics <- self$get_diagnostics(rp)
            }
            if (is.null(checks_made)) {
                checks_made <- list()
            }
            #Collecting Diagnostic Information
            # if (is.null(rpi$diagnostics)) {
            #     named_vars <-
            # }
            # named_vars <- rpi$diagnostics$diagnostic_information$names

            # summ <- posterior::summarize_draws(rp)
            # summ <- summ[which(summ[, 1] == named_vars), ]
            # rhat <- summ$rhat
            # ess_bulk <- summ$ess_bulk
            # ess_tail <- summ$ess_tail
            if (is.null(rpi$diagnostics)) {
                if (inherits(rp, "pdb_posterior_reference_draws")) {
                    diagnostics <- posterior::as_draws_df(rstan)
                }
                num_divergent <- rpi$diagnostics$divergent_transitions
            }
            if (sum(num_divergent) > 0) {
                stop(
                    "There were divergent transitions during sampling.",
                    call. = FALSE
                )
            }
            if (self$check_missing("ndraws_is_10k", checks_made)) {
                ndraws <- rpi$diagnostics$ndraws
                checks_made$ndraws_is_10k <- ndraws == 10000
            }
            if (self$check_missing("nchains_is_gte_4", checks_made)) {
                nchains <- rpi$diagnostics$nchains
                checks_made$nchains_is_gte_4 <- nchains >= 4
            }
            if (self$check_missing("r_hat_below_1_01", checks_made)) {
                rhat <- rpi$diagnostics$r_hat
                curr_len <- length(rhat)
                rhat_narm <- rhat |> na.omit()
                if (curr_len != length(rhat_narm)) {
                    message(
                        paste0(
                            "Some variables had undefined Rhat. ",
                            "This may be expected for constant or ",
                            "structurally constrained quantities. ",
                            "Verify that the corresponding draws ",
                            "are finite and legitimately constant "
                        )
                    )
                }
                checks_made$r_hat_below_1_01 <- self$check_rhat(rp, rhat)
            }

            if (self$check_missing("efmi_above_0_2", checks_made)) {
                efmi <- rpi$diagnostics$expected_fraction_of_missing_information
                checks_made$efmi_above_0_2 <- !anyNA(efmi) &&
                    all(is.finite(efmi)) &&
                    all(efmi > 0.2)
            }
            if (
                self$check_missing("abs_mean_lag1_ac_below_0_05", checks_made)
            ) {
                if (is.null(rpi$diagnostics$mean_lag1_ac)) {
                    mean_lag1_ac <- self$compute_mean_lag1_ac(
                        posterior::as_draws_array(rp)
                    )
                    rpi$diagnostics$mean_lag1_ac <- mean_lag1_ac
                    write_mean_ac <- TRUE
                } else {
                    mean_lag1_ac <- rpi$diagnostics$mean_lag1_ac
                }
                checks_made$abs_mean_lag1_ac_below_0_05 <- all(
                    mean_lag1_ac < 0.05,
                    na.rm = TRUE
                )
            }
            required_checks <- c(
                "ndraws_is_10k",
                "nchains_is_gte_4",
                "r_hat_below_1_01",
                "efmi_above_0_2",
                "abs_mean_lag1_ac_below_0_05"
            )
            checks_passed <- all(vapply(required_checks, function(x) {
                if (is.na(checks_made[[x]])) {
                    message(paste0(x, " could not be evaluated. "))
                    return(FALSE)
                } else if (!checks_made[[x]]) {
                    message(paste0(x, " did not pass."))
                    return(FALSE)
                }
                TRUE
            }, logical(1)))
            if (!checks_passed) {
                stop(call. = FALSE)
            }
            if (write_mean_ac && exists("mean_lag1_ac")) {
                info(rp)$diagnostics$mean_lag1_ac <-
                    rpi$diagnostics$mean_lag1_ac
            }
            info(rp)$checks_made <- checks_made
            invisible(rp)
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
            abs(sapply(var_names, function(name) {
                posterior::autocorrelation(
                    posterior::extract_variable(x, name)
                )[2]
            }))
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
                local_posterior_path,
                posterior_name
            ))
            jsonlite::read_json(posterior_path, simplifyVector = TRUE)
        },
        get_posterior_dims = function(po) {
            self$get_posterior_json(po)$dimensions
        },
        write_reference_draws = function(
            rp,
            recheck = TRUE,
            overwrite = FALSE,
            verify = TRUE
        ) {
            if (!inherits(rp, "pdb_reference_posterior_draws")) {
                stop(
                    "Reference draws must be of class `pdb_reference_posterior_draws`.",
                    call. = FALSE
                )
            }

            checkmate::assert_flag(recheck)
            checkmate::assert_flag(overwrite)
            checkmate::assert_flag(verify)

            if (recheck) {
                rp <- self$check_draws(rp)
            }

            rpi <- posteriordb::info(rp)
            required_checks <- c(
                "ndraws_is_10k",
                "nchains_is_gte_4",
                "r_hat_below_1_01",
                "abs_mean_lag1_ac_below_0_05",
                "efmi_above_0_2"
            )
            missing_checks <- setdiff(
                required_checks,
                names(rpi$checks_made)
            )
            if (length(missing_checks) > 0L) {
                stop(
                    paste0(
                        "Missing required reference-draw checks: ",
                        paste(missing_checks, collapse = ", ")
                    ),
                    call. = FALSE
                )
            }
            failed_checks <- required_checks[
                !vapply(
                    rpi$checks_made[required_checks],
                    isTRUE,
                    logical(1)
                )
            ]
            if (length(failed_checks) > 0L) {
                stop(
                    paste0(
                        "Reference draws failed: ",
                        paste(failed_checks, collapse = ", ")
                    ),
                    call. = FALSE
                )
            }
            divergences <- rpi$diagnostics$divergent_transitions
            if (
                anyNA(divergences) ||
                    any(!is.finite(divergences)) ||
                    sum(divergences) > 0
            ) {
                stop(
                    "Reference draws have invalid or divergent transitions.",
                    call. = FALSE
                )
            }

            pdb_root <- private$pdb$pdb_local_endpoint
            if (
                is.null(pdb_root) ||
                    length(pdb_root) != 1L ||
                    !dir.exists(pdb_root)
            ) {
                stop(
                    "Could not resolve the local PosteriorDB root.",
                    call. = FALSE
                )
            }

            info_dir <- file.path(
                pdb_root,
                "reference_posteriors",
                "draws",
                "info"
            )
            draws_dir <- file.path(
                pdb_root,
                "reference_posteriors",
                "draws",
                "draws"
            )
            dir.create(info_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(draws_dir, recursive = TRUE, showWarnings = FALSE)

            info_path <- file.path(
                info_dir,
                paste0(rpi$name, ".info.json")
            )
            draws_zip_path <- file.path(
                draws_dir,
                paste0(rpi$name, ".json.zip")
            )
            output_paths <- c(info_path, draws_zip_path)

            existing <- file.exists(output_paths)
            if (any(existing) && !overwrite) {
                stop(
                    paste0(
                        "Output already exists: ",
                        paste(output_paths[existing], collapse = ", ")
                    ),
                    call. = FALSE
                )
            }

            temp_info <- tempfile(
                pattern = paste0(rpi$name, "-"),
                tmpdir = info_dir,
                fileext = ".info.json"
            )
            temp_draws_dir <- tempfile(
                pattern = paste0(rpi$name, "-"),
                tmpdir = draws_dir
            )
            dir.create(temp_draws_dir)
            temp_draws <- file.path(
                temp_draws_dir,
                paste0(rpi$name, ".json")
            )
            temp_zip <- tempfile(
                pattern = paste0(rpi$name, "-"),
                tmpdir = draws_dir,
                fileext = ".json.zip"
            )
            on.exit(
                {
                    unlink(temp_info, force = TRUE)
                    unlink(temp_zip, force = TRUE)
                    unlink(temp_draws_dir, recursive = TRUE, force = TRUE)
                },
                add = TRUE
            )

            info_for_json <- rpi
            class(info_for_json) <- unique(c(class(info_for_json), "list"))
            info_json <- jsonlite::toJSON(
                info_for_json,
                pretty = TRUE,
                auto_unbox = TRUE,
                null = "null",
                digits = NA,
                encoding = "UTF-8"
            )
            draws_json <- jsonlite::toJSON(
                rp,
                pretty = TRUE,
                auto_unbox = TRUE,
                null = "null",
                digits = NA,
                encoding = "UTF-8"
            )
            Encoding(info_json) <- "UTF-8"
            Encoding(draws_json) <- "UTF-8"
            writeLines(info_json, temp_info, useBytes = TRUE)
            writeLines(draws_json, temp_draws, useBytes = TRUE)

            zip_status <- utils::zip(
                zipfile = temp_zip,
                files = temp_draws,
                flags = "-jq"
            )
            if (!identical(zip_status, 0L) || !file.exists(temp_zip)) {
                stop(
                    "Failed to create the reference-draw ZIP archive.",
                    call. = FALSE
                )
            }
            unlink(temp_draws, force = TRUE)

            if (!file.copy(temp_info, info_path, overwrite = overwrite)) {
                stop(
                    "Failed to write the reference-draw info JSON.",
                    call. = FALSE
                )
            }
            if (!file.copy(temp_zip, draws_zip_path, overwrite = overwrite)) {
                stop(
                    "Failed to write the reference-draw ZIP archive.",
                    call. = FALSE
                )
            }

            if (verify) {
                parsed_info <- jsonlite::fromJSON(
                    info_path,
                    simplifyVector = FALSE
                )
                if (!identical(parsed_info$name, rpi$name)) {
                    stop(
                        "Written reference-draw info failed verification.",
                        call. = FALSE
                    )
                }
                zip_listing <- utils::unzip(draws_zip_path, list = TRUE)
                expected_member <- paste0(rpi$name, ".json")
                if (
                    nrow(zip_listing) != 1L ||
                        basename(zip_listing$Name[1]) != expected_member
                ) {
                    stop(
                        "Written reference-draw archive failed verification.",
                        call. = FALSE
                    )
                }
            }

            invisible(rp)
        },
        is_constant = function(z, tolerance = sqrt(.Machine$double.eps)) {
            if (!all(is.finite(z))) {
                return(FALSE)
            }

            scale <- max(1, max(abs(z)))
            diff(range(z)) <= tolerance * scale
        },
        add_bibtex_entry = function(bibtex_str) {
            write(
                paste0("\n\n", bibtex_str, "\n"),
                file = private$get_reference_path(),
                append = TRUE
            )
            self$refresh()
        },
        write_reference_draw_info = function(
            info,
            overwrite = FALSE,
            local_refdraw_info = "/posterior_database/reference_posteriors/draws/info/"
        ) {
            if (inherits(info, "list")) {
                info <- posteriordb::as.pdb_reference_posterior_info(info)
            }
            if (!inherits(info, "pdb_reference_posterior_info")) {
                stop("Reference draw info needs to be posterior")
            }
            posteriordb::write_pdb(info, self$get_pdb(), overwrite = overwrite)
        },
        write_draws = function(draws) {
            checkmate::assertTRUE(inherits(draws, "draws"))
            checkmate::assertTRUE(inherits(
                draws,
                "pdb_reference_posterior_draws"
            ))
            posteriordb::write_pdb()
        },
        check_missing = function(name, checks) {
            value <- checks[[name]]
            is.null(value) || length(value) != 1L
        },
        add_bibtex_file = function(path) {},
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
            normalizePath(file.path(self$get_pdb_path(), local_model_path))
        },
        get_data_path = function(
            local_data_path = "/posterior_database/data/data/"
        ) {
            normalizePath(file.path(self$get_pdb_path(), local_data_path))
        },
        get_posterior_modeldata_files = function(posterior_name) {
            data_model_name <- strsplit(posterior_name, split = "-")[[1]]
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
        get_posterior_modeldata = function(posterior_name) {
            md_file <- self$get_posterior_modeldata_files(posterior_name)
            files_exist <- sapply(dm_paths, file.exists)
            if (!all(files_exist)) {
                stop(
                    paste0(
                        "The data or model file does not exist.",
                        "Please ensure you specified the correct posterior name."
                    ),
                    call. = FALSE
                )
            }
            posterior_model <- rstan::stan_model(
                file = dm_paths$model_file,
                model_name = posterior_name
            )
            posterior_model_code <- paste(
                readLines(md_file$data_file),
                collapse = "\n"
            )
            posterior_data <- unzip(self$copy_to_tempdir(
                dm_paths$data_file,
                unzip = TRUE,
                return_data = TRUE
            ))
            list(
                data = posterior_data,
                stan_model = posterior_model,
                model_code = posterior_model_code,
                stan_file = md_file$model_file
            )
        },
        # get_reference_posterior_info = function(
        #     po,
        #     local_posterior_info_path = "/posterior_database/reference_posteriors/draws/info/"
        # ) {
        #     rpi_path <- normalizePath(file.path(self$get_pdb_path(), )
        # },
        get_rpi_path = function(
            rpi_name,
            local_rpi_path = "/posterior_database/reference_posteriors/draws/info/"
        ) {
            rpi_name <- paste0(rpi_name, ".json")
            normalizePath(file.path(
                self$get_pdb_path(),
                local_rpi_path,
                rpi_name
            ))
        },
        get_rp_path = function(
            rp_name,
            local_rp_path = "/posterior_database/reference_posteriors/draws/draws/"
        ) {
            rp_name <- paste0(rp_name, ".json")
            normalizePath(file.path(
                self$get_pdb_path(),
                local_rp_path,
                rp_name
            ))
        },
        copy_to_tempdir = function(
            file_path,
            return_obj = TRUE,
            overwrite = TRUE
        ) {
            if (grepl(".zip", file_path)) {
                unzip <- TRUE
            }
            td <- normalizePath(base::tempdir())
            copied_path <- normalizePath(file.path(td, basename(file_path)))
            file.copy(
                from = file_path,
                to = td,
                overwrite = overwrite
            )
            if (return_obj && unzip) {
                jsonlite::read_json(unzip(copied_path), simplifyVector = TRUE)
            } else if (return_obj) {
                jsonlite::read_json(copied_path, simplifyVector = TRUE)
            } else {
                copied_path
            }
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
