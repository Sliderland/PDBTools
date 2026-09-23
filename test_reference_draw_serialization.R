# Run with: Rscript test_reference_draw_serialization.R
# No sampling or changes to an existing PosteriorDB repository are required.
source("PDBEntryBuilder_v2.R")

local({
  root <- tempfile("pdbtools-reference-serialization-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  posterior_name <- "toy-reference"
  variables <- c("alpha", "beta")
  values <- array(
    seq_len(4L * 3L * length(variables)) / 7,
    dim = c(4L, 3L, length(variables)),
    dimnames = list(
      iteration = as.character(seq_len(4L)),
      chain = as.character(seq_len(3L)),
      variable = variables
    )
  )
  draws <- posterior::as_draws_array(values)
  fit <- structure(list(), model_name = posterior_name)
  info_path <- file.path(root, paste0(posterior_name, ".info.json"))
  archive_path <- file.path(root, paste0(posterior_name, ".json.zip"))

  TestBuilder <- R6::R6Class(
    "ReferenceSerializationTestBuilder",
    inherit = PDBEntryBuilder,
    public = list(
      initialize = function() invisible(self),
      get_fit_model_name = function(fit) attr(fit, "model_name"),
      get_fit_draws = function(fit, format = "draws_array") draws,
      get_posterior_dims = function(po) setNames(rep(1L, 2L), variables),
      get_rpi_path = function(rpi_name, ...) info_path,
      get_rp_path = function(rp_name, ...) archive_path
    )
  )
  builder <- TestBuilder$new()

  jsonlite::write_json(
    list(
      name = posterior_name,
      diagnostics = list(
        diagnostic_information = list(names = variables)
      )
    ),
    info_path,
    auto_unbox = TRUE
  )

  write_archive <- function(object) {
    stage <- tempfile("reference-archive-")
    dir.create(stage)
    on.exit(unlink(stage, recursive = TRUE), add = TRUE)
    json_path <- file.path(stage, paste0(posterior_name, ".json"))
    jsonlite::write_json(object, json_path, digits = NA, null = "null")
    zip_path <- file.path(stage, "reference.json.zip")
    stopifnot(identical(
      utils::zip(zip_path, json_path, flags = "-jq"),
      0L
    ))
    stopifnot(file.copy(zip_path, archive_path, overwrite = TRUE))
  }

  builder$write_rpd_from_stan_fit(fit, verify = TRUE)
  written <- builder$read_reference_files(posterior_name)$draws
  stopifnot(
    identical(
      c(length(written), length(written[[1L]]),
        length(written[[1L]][[1L]])),
      c(3L, 2L, 4L)
    ),
    isTRUE(builder$verify_reference_files(fit))
  )

  # Existing array-oriented archives must remain readable and verifiable.
  write_archive(draws)
  stopifnot(isTRUE(builder$verify_reference_files(fit)))

  # A changed value must fail even when archive dimensions still match.
  damaged <- posterior::as_draws_list(draws)
  damaged[[1L]][[1L]][1L] <- -999
  write_archive(damaged)
  error <- tryCatch(
    builder$verify_reference_files(fit),
    error = identity
  )
  stopifnot(
    inherits(error, "error"),
    grepl("values changed", conditionMessage(error), fixed = TRUE)
  )

  # A truncated iteration vector must not be mistaken for a valid archive.
  truncated <- posterior::as_draws_list(draws)
  truncated[[1L]][[1L]] <- truncated[[1L]][[1L]][-1L]
  write_archive(truncated)
  error <- tryCatch(builder$verify_reference_files(fit), error = identity)
  stopifnot(
    inherits(error, "error"),
    grepl("not rectangular", conditionMessage(error), fixed = TRUE)
  )

  # The metadata supplies names, so its order must agree with the draws.
  write_archive(posterior::as_draws_list(draws))
  jsonlite::write_json(
    list(
      name = posterior_name,
      diagnostics = list(
        diagnostic_information = list(names = rev(variables))
      )
    ),
    info_path,
    auto_unbox = TRUE
  )
  error <- tryCatch(builder$verify_reference_files(fit), error = identity)
  stopifnot(
    inherits(error, "error"),
    grepl("variable names or order", conditionMessage(error), fixed = TRUE)
  )

  message("Reference-draw serialization regression test passed.")
})
