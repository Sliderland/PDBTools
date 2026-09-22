#!/usr/bin/env Rscript

# Check and repair posterior metadata entries where
# reference_posterior_name is null even though matching reference-draw files
# exist. A repair is made only when both the reference info JSON and draw
# archive exist and the reference info name agrees with the posterior name.

suppressPackageStartupMessages(library(jsonlite))

# Configuration -------------------------------------------------------------
database_roots <- c(
    "posterior_database",
    "../posteriordb/posterior_database"
)
database_root <- database_roots[dir.exists(database_roots)][1]
if (is.na(database_root)) {
    stop("Could not find the posterior_database directory.")
}

posterior_dir <- file.path(database_root, "posteriors")
draw_dir <- file.path(database_root, "reference_posteriors", "draws", "draws")
reference_info_dir <- file.path(
    database_root,
    "reference_posteriors",
    "draws",
    "info"
)

# Set to TRUE to write a CSV report; FALSE prints the matching records.
save_output <- FALSE
output_file <- "null_reference_name_with_draws.csv"

# Set to TRUE to repair safe matches automatically.
auto_fix <- TRUE

# Set to TRUE to search subdirectories below posterior_dir.
recursive <- TRUE

# Helpers -------------------------------------------------------------------
read_json_safely <- function(path) {
    tryCatch(
        jsonlite::fromJSON(path, simplifyVector = FALSE),
        error = function(e) NULL
    )
}

write_json_atomically <- function(object, path) {
    directory <- dirname(path)
    temporary <- tempfile(
        pattern = ".reference-link-",
        tmpdir = directory,
        fileext = ".json"
    )
    backup <- tempfile(
        pattern = ".reference-link-backup-",
        tmpdir = directory,
        fileext = ".json"
    )
    on.exit(
        unlink(c(temporary, backup), force = TRUE),
        add = TRUE
    )

    jsonlite::write_json(
        object,
        temporary,
        pretty = TRUE,
        auto_unbox = TRUE,
        null = "null",
        digits = NA
    )

    had_original <- file.exists(path)
    if (had_original && !file.rename(path, backup)) {
        stop("Could not stage the existing posterior JSON: ", path)
    }
    if (!file.rename(temporary, path)) {
        if (had_original) {
            file.rename(backup, path)
        }
        stop("Could not install the repaired posterior JSON: ", path)
    }
    if (had_original) {
        unlink(backup, force = TRUE)
    }
    invisible(path)
}

posterior_files <- list.files(
    posterior_dir,
    pattern = "\\.json$",
    full.names = TRUE,
    recursive = recursive,
    ignore.case = TRUE
)

results <- lapply(posterior_files, function(path) {
    posterior <- read_json_safely(path)

    if (!is.list(posterior) ||
        !"reference_posterior_name" %in% names(posterior)) {
        return(NULL)
    }

    posterior_name <- posterior[["name"]]
    if (is.null(posterior_name) || length(posterior_name) != 1L) {
        posterior_name <- sub("\\.json$", "", basename(path), ignore.case = TRUE)
    }
    if (length(posterior_name) != 1L ||
        !is.character(posterior_name) ||
        is.na(posterior_name) ||
        !nzchar(posterior_name)) {
        return(NULL)
    }
    posterior_name <- as.character(posterior_name)

    explicit_reference_name <- posterior[["reference_posterior_name"]]
    is_null_link <- is.null(explicit_reference_name)
    if (!is_null_link &&
        (length(explicit_reference_name) != 1L ||
         !is.character(explicit_reference_name) ||
         is.na(explicit_reference_name) ||
         !nzchar(explicit_reference_name))) {
        return(data.frame(
            posterior_name = posterior_name,
            reference_posterior_name = NA_character_,
            expected_reference_name = NA_character_,
            posterior_file = path,
            reference_info_file = NA_character_,
            reference_draw_file = NA_character_,
            reference_info_exists = FALSE,
            reference_draw_exists = FALSE,
            reference_info_name = NA_character_,
            info_name_matches = FALSE,
            issue_type = "invalid_reference_name",
            auto_fix = auto_fix,
            fixed = FALSE,
            fix_status = "not_auto_fixed",
            fix_error = NA_character_,
            stringsAsFactors = FALSE
        ))
    }

    reference_name <- if (is_null_link) {
        posterior_name
    } else {
        as.character(explicit_reference_name)
    }
    draw_path <- file.path(draw_dir, paste0(reference_name, ".json.zip"))
    info_path <- file.path(
        reference_info_dir,
        paste0(reference_name, ".info.json")
    )
    draw_exists <- file.exists(draw_path)
    info_exists <- file.exists(info_path)
    reference_info <- if (info_exists) {
        read_json_safely(info_path)
    } else {
        NULL
    }
    reference_info_name <- if (
        is.list(reference_info) &&
            !is.null(reference_info$name) &&
            length(reference_info$name) == 1L
    ) {
        as.character(reference_info$name)
    } else {
        NA_character_
    }
    info_name_matches <- !is.na(reference_info_name) &&
        identical(reference_info_name, reference_name)
    complete_reference <- draw_exists && info_exists && info_name_matches

    # A null link is valid when no reference material exists. If reference
    # files do exist, report the mismatch and repair only a complete,
    # canonical pair. Explicit non-null links are never guessed or rewritten.
    if (is_null_link) {
        if (!draw_exists && !info_exists) {
            return(NULL)
        }
        issue_type <- if (complete_reference) {
            "null_link_with_reference_files"
        } else if (!draw_exists) {
            "missing_reference_draw"
        } else if (!info_exists) {
            "missing_reference_info"
        } else if (is.null(reference_info)) {
            "invalid_reference_info"
        } else {
            "reference_info_name_mismatch"
        }
    } else {
        if (complete_reference) {
            return(NULL)
        }
        issue_type <- if (!draw_exists) {
            "missing_reference_draw"
        } else if (!info_exists) {
            "missing_reference_info"
        } else if (is.null(reference_info)) {
            "invalid_reference_info"
        } else {
            "reference_info_name_mismatch"
        }
    }

    fixed <- FALSE
    fix_status <- if (is_null_link && complete_reference) "fixable" else "not_auto_fixed"
    fix_error <- NA_character_

    if (auto_fix && is_null_link && complete_reference) {
        tryCatch(
            {
                posterior[["reference_posterior_name"]] <- reference_name
                write_json_atomically(posterior, path)
                written <- read_json_safely(path)
                if (
                    !is.list(written) ||
                        !identical(
                            written[["reference_posterior_name"]],
                            reference_name
                        )
                ) {
                    stop("Round-trip verification failed.")
                }
                fixed <- TRUE
                fix_status <- "fixed"
            },
            error = function(e) {
                fix_status <<- "fix_failed"
                fix_error <<- conditionMessage(e)
            }
        )
    }

    data.frame(
        posterior_name = posterior_name,
        reference_posterior_name = if (is_null_link) NA_character_ else reference_name,
        expected_reference_name = reference_name,
        posterior_file = path,
        reference_info_file = info_path,
        reference_draw_file = draw_path,
        reference_info_exists = info_exists,
        reference_draw_exists = draw_exists,
        reference_info_name = reference_info_name,
        info_name_matches = info_name_matches,
        issue_type = issue_type,
        auto_fix = auto_fix,
        fixed = fixed,
        fix_status = fix_status,
        fix_error = fix_error,
        stringsAsFactors = FALSE
    )
})

results <- Filter(Negate(is.null), results)

if (length(results) == 0L) {
    results <- data.frame(
        posterior_name = character(),
        reference_posterior_name = character(),
        expected_reference_name = character(),
        posterior_file = character(),
        reference_info_file = character(),
        reference_draw_file = character(),
        reference_info_exists = logical(),
        reference_draw_exists = logical(),
        reference_info_name = character(),
        info_name_matches = logical(),
        issue_type = character(),
        auto_fix = logical(),
        fixed = logical(),
        fix_status = character(),
        fix_error = character(),
        stringsAsFactors = FALSE
    )
} else {
    results <- do.call(rbind, results)
}

# Every row is an inconsistency; consistent explicit links were omitted above.
matches <- results

if (save_output) {
    write.csv(matches, output_file, row.names = FALSE)
    message("Saved ", nrow(matches), " inconsistent posterior(s) to ", output_file)
} else if (nrow(matches) == 0L) {
    message("No inconsistent reference_posterior_name entries were found.")
} else {
    print(matches, row.names = FALSE)
}

invisible(matches)
