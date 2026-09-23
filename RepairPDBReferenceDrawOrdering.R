#!/usr/bin/env Rscript

# Audit and, optionally, repair reference-draw JSON archives written in the
# iteration x chain x variable layout.  PosteriorDB's standard JSON layout is
# chain x variable x iteration.
#
# Dry run:
#   Rscript RepairPDBReferenceDrawOrdering.R
#
# Audit one archive:
#   Rscript RepairPDBReferenceDrawOrdering.R \
#     --posterior-name=heaps_med10-statprior_var
#
# Apply repairs (creates a .bak-TIMESTAMP archive first):
#   Rscript RepairPDBReferenceDrawOrdering.R --apply
#
# Write a repaired, standalone copy instead of changing the source tree:
#   Rscript RepairPDBReferenceDrawOrdering.R \
#     --apply --output-pdb-path=~/Documents/posteriordb_repaired

suppressPackageStartupMessages({
  library(jsonlite)
  library(posterior)
})

parse_args <- function(args) {
  options <- list(
    pdb_path = path.expand("~/Documents/posteriordb"),
    output_pdb_path = NULL,
    posterior_name = NULL,
    apply = FALSE,
    backup = TRUE,
    report = NULL
  )

  for (arg in args) {
    if (identical(arg, "--apply")) {
      options$apply <- TRUE
    } else if (identical(arg, "--no-backup")) {
      options$backup <- FALSE
    } else if (grepl("^--pdb-path=", arg)) {
      options$pdb_path <- path.expand(sub("^--pdb-path=", "", arg))
    } else if (grepl("^--output-pdb-path=", arg)) {
      options$output_pdb_path <- path.expand(
        sub("^--output-pdb-path=", "", arg)
      )
    } else if (grepl("^--posterior-name=", arg)) {
      options$posterior_name <- sub("^--posterior-name=", "", arg)
    } else if (grepl("^--report=", arg)) {
      options$report <- path.expand(sub("^--report=", "", arg))
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  if (is.null(options$report)) {
    report_root <- if (is.null(options$output_pdb_path)) {
      options$pdb_path
    } else {
      options$output_pdb_path
    }
    options$report <- file.path(
      report_root,
      "reference_draw_order_audit.csv"
    )
  }

  options
}

copy_pdb_tree <- function(source_path, destination_path) {
  if (dir.exists(destination_path)) {
    existing <- list.files(
      destination_path,
      all.files = TRUE,
      no.. = TRUE,
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(existing) > 0L) {
      stop(
        "Output PDB path already exists and is not empty: ",
        destination_path,
        call. = FALSE
      )
    }
  } else if (!dir.create(destination_path, recursive = TRUE)) {
    stop("Could not create output PDB path: ", destination_path, call. = FALSE)
  }

  files <- list.files(
    source_path,
    all.files = TRUE,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  files <- files[!file.info(files)$isdir]
  for (source_file in files) {
    relative_file <- substring(
      source_file,
      nchar(normalizePath(source_path, mustWork = TRUE)) + 2L
    )
    destination_file <- file.path(destination_path, relative_file)
    dir.create(dirname(destination_file), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source_file, destination_file, overwrite = FALSE)) {
      stop("Could not copy: ", source_file, call. = FALSE)
    }
  }
  invisible(destination_path)
}

reference_root <- function(pdb_path) {
  file.path(
    pdb_path,
    "posterior_database",
    "reference_posteriors"
  )
}

draw_archive_root <- function(pdb_path) {
  file.path(reference_root(pdb_path), "draws", "draws")
}

draw_info_root <- function(pdb_path) {
  file.path(reference_root(pdb_path), "draws", "info")
}

summary_root <- function(pdb_path) {
  file.path(reference_root(pdb_path), "summary_statistics")
}

read_draw_archive <- function(archive_path) {
  members <- unzip(archive_path, list = TRUE)$Name
  json_member <- members[
    grepl("\\.json$", members, ignore.case = TRUE)
  ]

  if (length(json_member) != 1L) {
    stop(
      "Expected one JSON member in archive: ",
      archive_path,
      call. = FALSE
    )
  }

  extraction_dir <- tempfile("pdb_draw_archive_")
  dir.create(extraction_dir)
  on.exit(unlink(extraction_dir, recursive = TRUE), add = TRUE)
  extracted <- unzip(
    archive_path,
    files = json_member,
    exdir = extraction_dir,
    junkpaths = TRUE
  )

  jsonlite::read_json(extracted, simplifyVector = FALSE)
}

read_reference_metadata <- function(info_path) {
  info <- jsonlite::read_json(info_path, simplifyVector = FALSE)
  variable_names <- info$diagnostics$diagnostic_information$names
  number_of_chains <- as.integer(
    info$inference$method_arguments$chains
  )

  if (is.null(variable_names) || length(variable_names) == 0L) {
    stop("Missing variable names in: ", info_path, call. = FALSE)
  }
  if (length(number_of_chains) != 1L || is.na(number_of_chains)) {
    stop("Missing chain count in: ", info_path, call. = FALSE)
  }

  list(
    info = info,
    variable_names = unname(unlist(variable_names, use.names = FALSE)),
    number_of_chains = number_of_chains
  )
}

raw_dimensions <- function(raw_draws) {
  c(
    length(raw_draws),
    length(raw_draws[[1L]]),
    length(raw_draws[[1L]][[1L]])
  )
}

classify_layout <- function(raw_draws, metadata) {
  dimensions <- raw_dimensions(raw_draws)
  number_of_chains <- metadata$number_of_chains
  number_of_variables <- length(metadata$variable_names)

  standard <- dimensions[1L] == number_of_chains &&
    dimensions[2L] == number_of_variables
  array_order <- dimensions[2L] == number_of_chains &&
    dimensions[3L] == number_of_variables

  if (standard && !array_order) {
    return("chain_variable_iteration")
  }
  if (array_order && !standard) {
    return("iteration_chain_variable")
  }
  if (standard && array_order) {
    return("ambiguous_symmetric")
  }
  "unrecognized"
}

raw_to_draws_array <- function(raw_draws, metadata, layout) {
  dimensions <- raw_dimensions(raw_draws)
  values <- unlist(raw_draws, recursive = TRUE, use.names = FALSE)

  if (identical(layout, "chain_variable_iteration")) {
    values <- aperm(
      array(
        as.numeric(values),
        dim = c(dimensions[3L], dimensions[2L], dimensions[1L])
      ),
      perm = c(1L, 3L, 2L)
    )
  } else if (identical(layout, "iteration_chain_variable")) {
    values <- aperm(
      array(
        as.numeric(values),
        dim = dimensions[c(3L, 2L, 1L)]
      ),
      perm = c(3L, 2L, 1L)
    )
  } else {
    stop("Cannot normalize layout: ", layout, call. = FALSE)
  }

  dimnames(values) <- list(
    iteration = as.character(seq_len(dim(values)[1L])),
    chain = as.character(seq_len(dim(values)[2L])),
    variable = metadata$variable_names
  )
  posterior::as_draws_array(values)
}

write_standard_archive <- function(draws_array, archive_path) {
  draws_list <- posterior::as_draws_list(draws_array)
  temporary_dir <- tempfile("pdb_repaired_archive_")
  dir.create(temporary_dir)
  on.exit(unlink(temporary_dir, recursive = TRUE), add = TRUE)

  json_path <- file.path(
    temporary_dir,
    sub("\\.zip$", "", basename(archive_path))
  )
  # Keep the output archive outside the staging directory.  Some platform
  # zip implementations cannot create the archive in the same directory
  # while using absolute input paths and the -j (junk paths) flag.
  zip_path <- tempfile(
    pattern = "pdb_repaired_archive_",
    fileext = ".zip"
  )

  jsonlite::write_json(
    draws_list,
    json_path,
    digits = NA,
    null = "null"
  )
  status <- utils::zip(
    zipfile = zip_path,
    files = json_path,
    flags = "-jq"
  )
  if (!identical(status, 0L) || !file.exists(zip_path)) {
    stop("Could not create repaired archive.", call. = FALSE)
  }

  zip_path
}

audit_summary_file <- function(path) {
  summary <- jsonlite::read_json(path, simplifyVector = FALSE)
  names <- summary$names
  if (is.null(names)) {
    return(list(status = "missing_names", variables = NA_integer_))
  }

  names <- unlist(names, use.names = FALSE)
  fields <- setdiff(names(summary), "names")
  lengths_ok <- all(vapply(
    summary[fields],
    function(x) length(x) == length(names),
    logical(1)
  ))

  list(
    status = if (lengths_ok) "ok" else "length_mismatch",
    variables = length(names)
  )
}

audit_summaries <- function(pdb_path, posterior_name = NULL) {
  root <- summary_root(pdb_path)
  paths <- list.files(
    root,
    pattern = "\\.json$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (!is.null(posterior_name)) {
    paths <- paths[basename(paths) == paste0(posterior_name, ".json")]
  }
  if (length(paths) == 0L) {
    return(data.frame())
  }

  rows <- lapply(paths, function(path) {
    result <- tryCatch(
      audit_summary_file(path),
      error = function(error) {
        list(status = paste0("error: ", conditionMessage(error)),
             variables = NA_integer_)
      }
    )
    data.frame(
      file = path,
      status = result$status,
      variables = result$variables,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

audit_and_maybe_repair <- function(
    pdb_path,
    posterior_name = NULL,
    apply = FALSE,
    backup = TRUE
) {
  archive_paths <- list.files(
    draw_archive_root(pdb_path),
    pattern = "\\.json\\.zip$",
    full.names = TRUE
  )
  if (!is.null(posterior_name)) {
    archive_paths <- archive_paths[
      basename(archive_paths) == paste0(posterior_name, ".json.zip")
    ]
  }
  if (length(archive_paths) == 0L) {
    stop("No matching reference-draw archives were found.", call. = FALSE)
  }

  rows <- lapply(archive_paths, function(archive_path) {
    name <- sub("\\.json\\.zip$", "", basename(archive_path))
    info_path <- file.path(draw_info_root(pdb_path), paste0(name, ".info.json"))
    result <- data.frame(
      posterior_name = name,
      raw_dimensions = NA_character_,
      chains = NA_integer_,
      variables = NA_integer_,
      layout = NA_character_,
      action = "error",
      message = NA_character_,
      stringsAsFactors = FALSE
    )

    tryCatch({
      metadata <- read_reference_metadata(info_path)
      raw_draws <- read_draw_archive(archive_path)
      dimensions <- raw_dimensions(raw_draws)
      layout <- classify_layout(raw_draws, metadata)
      result$raw_dimensions <- paste(dimensions, collapse = " x ")
      result$chains <- metadata$number_of_chains
      result$variables <- length(metadata$variable_names)
      result$layout <- layout

      if (identical(layout, "chain_variable_iteration")) {
        result$action <- "already_standard"
        result$message <- "No rewrite needed."
      } else if (identical(layout, "iteration_chain_variable")) {
        draws_array <- raw_to_draws_array(raw_draws, metadata, layout)
        repaired_zip <- write_standard_archive(draws_array, archive_path)

        repaired_raw <- read_draw_archive(repaired_zip)
        repaired_layout <- classify_layout(repaired_raw, metadata)
        if (!identical(repaired_layout, "chain_variable_iteration")) {
          stop("Repaired archive did not validate as standard layout.")
        }

        if (!apply) {
          result$action <- "would_rewrite"
          result$message <- "Dry run; archive unchanged."
        } else {
          timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
          backup_path <- paste0(archive_path, ".bak-", timestamp)
          if (backup && !file.copy(archive_path, backup_path)) {
            stop("Could not create backup: ", backup_path)
          }
          if (!file.copy(repaired_zip, archive_path, overwrite = TRUE)) {
            stop("Could not replace archive: ", archive_path)
          }
          result$action <- "rewritten"
          result$message <- if (backup) {
            paste("Backup:", backup_path)
          } else {
            "Archive rewritten without backup."
          }
        }
        unlink(repaired_zip, force = TRUE)
      } else {
        result$action <- "manual_review"
        result$message <- paste("Layout was", layout)
      }
    }, error = function(error) {
      result$message <<- conditionMessage(error)
    })
    result
  })

  do.call(rbind, rows)
}

main <- function() {
  options <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!dir.exists(options$pdb_path)) {
    stop("PosteriorDB path does not exist: ", options$pdb_path, call. = FALSE)
  }
  if (!is.null(options$output_pdb_path) && !options$apply) {
    stop(
      "--output-pdb-path requires --apply; use a dry run without it first.",
      call. = FALSE
    )
  }

  working_pdb_path <- options$pdb_path
  if (!is.null(options$output_pdb_path)) {
    source_path <- normalizePath(options$pdb_path, mustWork = TRUE)
    output_parent <- dirname(options$output_pdb_path)
    dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
    output_path <- normalizePath(
      options$output_pdb_path,
      mustWork = FALSE
    )
    if (identical(source_path, output_path)) {
      stop("Output PDB path must differ from the source path.", call. = FALSE)
    }
    message("Copying source PosteriorDB tree to: ", output_path)
    copy_pdb_tree(source_path, output_path)
    working_pdb_path <- output_path
  }

  message(
    if (!is.null(options$output_pdb_path)) {
      "Applying repairs to the copied PosteriorDB tree."
    } else if (options$apply) {
      "Applying reference-draw repairs in place."
    } else {
      "Dry run: no reference-draw archives will be changed."
    }
  )
  draw_report <- audit_and_maybe_repair(
    pdb_path = working_pdb_path,
    posterior_name = options$posterior_name,
    apply = options$apply,
    backup = options$backup
  )
  summary_report <- audit_summaries(
    pdb_path = working_pdb_path,
    posterior_name = options$posterior_name
  )

  print(draw_report)
  if (nrow(summary_report) == 0L) {
    message("No matching summary-statistic files found.")
  } else {
    message("Summary-statistic audit:")
    print(as.data.frame(table(summary_report$status)))
  }

  report <- options$report
  dir.create(dirname(report), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(draw_report, report, row.names = FALSE)
  message("Draw audit written to: ", report)

  summary_report_path <- sub(
    "\\.csv$",
    "_summaries.csv",
    report,
    ignore.case = TRUE
  )
  if (!identical(summary_report_path, report)) {
    utils::write.csv(summary_report, summary_report_path, row.names = FALSE)
    message("Summary audit written to: ", summary_report_path)
  }
}

if (identical(environment(), globalenv())) {
  main()
}
