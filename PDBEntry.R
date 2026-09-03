library(R6)

.pdb_result <- function(status = "pending", object = NULL, error = NULL) {
    list(
        status = status,
        object = object,
        error = if (is.null(error)) {
            NULL
        } else if (inherits(error, "condition")) {
            conditionMessage(error)
        } else {
            as.character(error)
        }
    )
}
.pdb_name <- function(x, label = "name") {
    if (!is.character(x) || length(x) != 1L || is.na(x) || x == "") {
        stop("`", label, "` must be one non-empty string.", call. = FALSE)
    }
}
.pdb_set <- function(x, stage, status, object = NULL, error = NULL) {
    if (
        !status %in%
            c("pending", "valid", "invalid", "added", "skipped", "failed")
    ) {
        stop("Unknown status.", call. = FALSE)
    }
    x[[stage]] <- .pdb_result(status, object, error)
    x
}
.pdb_results <- function(stages) {
    setNames(lapply(stages, function(x) .pdb_result()), stages)
}

PDBDataEntry <- R6Class(
    "PDBDataEntry",
    public = list(
        name = NULL,
        raw_data = NULL,
        data_info = NULL,
        pdb_data = NULL,
        results = NULL,
        initialize = function(name, raw_data, data_info) {
            .pdb_name(name)
            if (!is.list(raw_data)) {
                stop("`raw_data` must be a named list.", call. = FALSE)
            }
            if (!is.list(data_info) && !inherits(data_info, "pdb_data_info")) {
                stop("Invalid `data_info`.", call. = FALSE)
            }
            n <- tryCatch(data_info$name, error = function(e) NULL)
            if (!is.null(n) && as.character(n) != name) {
                stop("Name differs from `data_info$name`.", call. = FALSE)
            }
            if (is.list(data_info) && is.null(data_info$name)) {
                data_info$name <- name
            }
            self$name <- name
            self$raw_data <- raw_data
            self$data_info <- data_info
            self$reset()
            invisible(self)
        },
        reset = function() {
            self$pdb_data <- NULL
            self$results <- .pdb_results(c("validation", "creation", "writing"))
            invisible(self)
        },
        set_result = function(stage, status, object = NULL, error = NULL) {
            self$results <- .pdb_set(self$results, stage, status, object, error)
            if (!is.null(object) && stage %in% c("creation", "writing")) {
                self$pdb_data <- object
            }
            invisible(self)
        },
        as_list = function() {
            list(
                name = self$name,
                raw_data = self$raw_data,
                data_info = self$data_info,
                pdb_data = self$pdb_data,
                results = self$results
            )
        }
    )
)

PDBModelEntry <- R6Class(
    "PDBModelEntry",
    public = list(
        name = NULL,
        stan_file = NULL,
        model_info = NULL,
        pdb_model_code = NULL,
        results = NULL,
        initialize = function(name, stan_file, model_info) {
            .pdb_name(name)
            if (
                !is.character(stan_file) ||
                    length(stan_file) != 1L ||
                    is.na(stan_file) ||
                    stan_file == ""
            ) {
                stop("Invalid `stan_file`.", call. = FALSE)
            }
            if (
                !is.list(model_info) && !inherits(model_info, "pdb_model_info")
            ) {
                stop("Invalid `model_info`.", call. = FALSE)
            }
            n <- tryCatch(model_info$name, error = function(e) NULL)
            if (!is.null(n) && as.character(n) != name) {
                stop("Name differs from `model_info$name`.", call. = FALSE)
            }
            if (is.list(model_info) && is.null(model_info$name)) {
                model_info$name <- name
            }
            self$name <- name
            self$stan_file <- stan_file
            self$model_info <- model_info
            self$reset()
            invisible(self)
        },
        reset = function() {
            self$pdb_model_code <- NULL
            self$results <- .pdb_results(c("validation", "creation", "writing"))
            invisible(self)
        },
        set_result = function(stage, status, object = NULL, error = NULL) {
            self$results <- .pdb_set(self$results, stage, status, object, error)
            if (!is.null(object) && stage %in% c("creation", "writing")) {
                self$pdb_model_code <- object
            }
            invisible(self)
        },
        as_list = function() {
            list(
                name = self$name,
                stan_file = self$stan_file,
                model_info = self$model_info,
                pdb_model_code = self$pdb_model_code,
                results = self$results
            )
        }
    )
)

PDBPosteriorEntry <- R6Class(
    "PDBPosteriorEntry",
    public = list(
        name = NULL,
        data_entry = NULL,
        model_entry = NULL,
        posterior_spec = NULL,
        pdb_posterior = NULL,
        reference_sampling = NULL,
        reference_info = NULL,
        reference_draws = NULL,
        results = NULL,
        initialize = function(
            name,
            data_entry,
            model_entry,
            posterior_spec = list(),
            reference_sampling = NULL
        ) {
            .pdb_name(name)
            if (!inherits(data_entry, "PDBDataEntry")) {
                stop("Invalid `data_entry`.", call. = FALSE)
            }
            if (!inherits(model_entry, "PDBModelEntry")) {
                stop("Invalid `model_entry`.", call. = FALSE)
            }
            if (!is.list(posterior_spec)) {
                stop("`posterior_spec` must be a list.", call. = FALSE)
            }
            if (!is.null(reference_sampling) && !is.list(reference_sampling)) {
                stop("Invalid `reference_sampling`.", call. = FALSE)
            }
            self$name <- name
            self$data_entry <- data_entry
            self$model_entry <- model_entry
            self$posterior_spec <- posterior_spec
            self$reference_sampling <- reference_sampling
            self$reset()
            invisible(self)
        },
        reset = function() {
            self$pdb_posterior <- NULL
            self$reference_info <- NULL
            self$reference_draws <- NULL
            self$results <- .pdb_results(c(
                "validation",
                "creation",
                "writing",
                "reference"
            ))
            invisible(self)
        },
        set_result = function(stage, status, object = NULL, error = NULL) {
            self$results <- .pdb_set(self$results, stage, status, object, error)
            if (!is.null(object) && stage %in% c("creation", "writing")) {
                self$pdb_posterior <- object
            }
            if (!is.null(object) && stage == "reference") {
                self$reference_draws <- object
            }
            invisible(self)
        },
        set_reference_info = function(x) {
            self$reference_info <- x
            invisible(self)
        },
        as_list = function() {
            list(
                name = self$name,
                data_name = self$data_entry$name,
                model_name = self$model_entry$name,
                data_entry = self$data_entry,
                model_entry = self$model_entry,
                posterior_spec = self$posterior_spec,
                pdb_posterior = self$pdb_posterior,
                reference_sampling = self$reference_sampling,
                reference_info = self$reference_info,
                reference_draws = self$reference_draws,
                results = self$results
            )
        }
    )
)

PDBEntry <- PDBPosteriorEntry
