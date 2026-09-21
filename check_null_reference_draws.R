#!/usr/bin/env Rscript

# Check for posterior metadata entries where reference_posterior_name is null
# even though a matching reference-draw archive exists.

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

# Set to TRUE to write a CSV report; FALSE prints matching files instead.
save_output <- TRUE
output_file <- "null_reference_name_with_draws.csv"

# Set to TRUE to search subdirectories below posterior_dir.
recursive <- TRUE

# Helpers -------------------------------------------------------------------
read_posterior <- function(path) {
    tryCatch(
        jsonlite::fromJSON(path, simplifyVector = FALSE),
        error = function(e) NULL
    )
}

posterior_files <- list.files(
    posterior_dir,
    pattern = "\\.json$",
    full.names = TRUE,
    recursive = recursive,
    ignore.case = TRUE
)

results <- lapply(posterior_files, function(path) {
    posterior <- read_posterior(path)

    if (!is.list(posterior) ||
        !"reference_posterior_name" %in% names(posterior) ||
        !is.null(posterior[["reference_posterior_name"]])) {
        return(NULL)
    }

    posterior_name <- posterior[["name"]]
    if (is.null(posterior_name) || length(posterior_name) != 1L) {
        posterior_name <- sub("\\.json$", "", basename(path), ignore.case = TRUE)
    }

    draw_path <- file.path(draw_dir, paste0(posterior_name, ".json.zip"))

    data.frame(
        posterior_name = posterior_name,
        posterior_file = path,
        reference_draw_file = draw_path,
        reference_draw_exists = file.exists(draw_path),
        stringsAsFactors = FALSE
    )
})

results <- Filter(Negate(is.null), results)

if (length(results) == 0L) {
    results <- data.frame(
        posterior_name = character(),
        posterior_file = character(),
        reference_draw_file = character(),
        reference_draw_exists = logical(),
        stringsAsFactors = FALSE
    )
} else {
    results <- do.call(rbind, results)
}

# Keep only the cases that appear inconsistent: null metadata plus draws.
matches <- results[results$reference_draw_exists, , drop = FALSE]

if (save_output) {
    write.csv(matches, output_file, row.names = FALSE)
    message("Saved ", nrow(matches), " matching posterior(s) to ", output_file)
} else if (nrow(matches) == 0L) {
    message("No null reference_posterior_name entries have matching draws.")
} else {
    writeLines(matches$posterior_file)
}

invisible(matches)
