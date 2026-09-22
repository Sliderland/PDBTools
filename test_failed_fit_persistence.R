# Exercise per-workflow fit persistence without running a Stan sampler.
source("PDBEntryBuilder_v2.R")

local({
  diagnostic_dir <- tempfile("failed_fit_test_")
  on.exit(unlink(diagnostic_dir, recursive = TRUE), add = TRUE)

  stopifnot(requireNamespace("rstan", quietly = TRUE))
  fit <- methods::new("stanfit")
  test_mode <- "error"
  write_calls <- character()
  TestBuilder <- R6::R6Class(
    "TestBuilder",
    inherit = PDBEntryBuilder,
    public = list(
      compute_reference_draws = function(on_sampled, ...) {
        sampled_fit <- if (identical(test_mode, "unsupported")) list() else fit
        on_sampled(sampled_fit)
        if (identical(test_mode, "interrupt")) {
          stop(structure(
            list(message = "simulated interrupt", call = NULL),
            class = c("interrupt", "condition")
          ))
        }
        sampled_fit
      },
      get_diagnostics = function(...) {
        if (test_mode %in% c("error", "unsupported")) {
          stop("simulated diagnostic failure")
        }
        list(divergent_transitions = if (identical(test_mode, "divergent")) 1 else 0)
      },
      get_checks_from_stanfit = function(...) {
        list(
          ndraws_is_10k = TRUE,
          nchains_is_gte_4 = TRUE,
          r_hat_below_1_01 = TRUE,
          efmi_above_0_2 = TRUE,
          abs_mean_lag1_ac_below_0_05 = TRUE
        )
      },
      check_draws_from_stanfit = function(stan_fit) stan_fit,
      write_rpi_from_stan_fit = function(...) {
        write_calls <<- c(write_calls, "info")
        "test-info.json"
      },
      write_rpd_from_stan_fit = function(...) {
        write_calls <<- c(write_calls, "draws")
        "test-draws.csv"
      },
      write_summary_statistics_from_stan_fit = function(...) {
        write_calls <<- c(write_calls, "summary")
        "test-summary.json"
      },
      verify_reference_files = function(...) {
        write_calls <<- c(write_calls, "verify")
        invisible(TRUE)
      },
      link_reference_posterior_from_stan_fit = function(...) {
        write_calls <<- c(write_calls, "link")
        invisible(TRUE)
      }
    )
  )
  builder <- TestBuilder$new(
    path.expand("~/Documents/posteriordb"),
    detect_cores = FALSE,
    n_cores = 1L
  )
  entry_spec <- list(
    posterior_name = "test-failed-fit",
    posterior = list(data_name = "test-data", model_name = "test-model"),
    reference = list(sampling_args = list())
  )
  run_case <- function(write = FALSE) {
    builder$run_workflows(
      entries = list(test = entry_spec),
      register = FALSE,
      write = write,
      save_failed_fits = TRUE,
      failed_fit_dir = diagnostic_dir
    )[[1L]]
  }

  failed <- run_case()
  stopifnot(
    identical(failed$status, "error"),
    identical(failed$error, "simulated diagnostic failure"),
    file.exists(failed$failed_fit_path),
    inherits(readRDS(failed$failed_fit_path), "stanfit")
  )
  metadata_path <- sub("\\.rds$", ".diagnostics.rds", failed$failed_fit_path)
  stopifnot(
    file.exists(metadata_path),
    identical(readRDS(metadata_path)$reason$status, "error")
  )

  test_mode <- "divergent"
  bad_checks <- run_case(write = TRUE)
  stopifnot(
    identical(bad_checks$status, "failed_checks"),
    identical(bad_checks$failed_required_checks, "divergent_transitions"),
    !isTRUE(bad_checks$written),
    file.exists(bad_checks$failed_fit_path),
    identical(write_calls, character())
  )
  bad_metadata <- sub("\\.rds$", ".diagnostics.rds", bad_checks$failed_fit_path)
  stopifnot(identical(readRDS(bad_metadata)$reason$status, "failed_checks"))

  test_mode <- "pass"
  passing <- run_case(write = TRUE)
  stopifnot(
    identical(passing$status, "completed"),
    isTRUE(passing$written),
    inherits(passing$fit, "stanfit"),
    is.null(passing$failed_fit_path),
    length(list.files(diagnostic_dir, pattern = "\\.rds$")) == 4L,
    identical(write_calls, c("info", "draws", "summary", "verify", "link"))
  )

  test_mode <- "interrupt"
  interrupted <- tryCatch(run_case(), interrupt = function(e) e)
  stopifnot(inherits(interrupted, "interrupt"))
  metadata <- list.files(
    diagnostic_dir, pattern = "\\.diagnostics\\.rds$", full.names = TRUE
  )
  stopifnot(
    length(metadata) == 3L,
    any(vapply(metadata, function(path) {
      identical(readRDS(path)$reason$status, "interrupted")
    }, logical(1)))
  )

  test_mode <- "unsupported"
  unsupported <- tryCatch(run_case(), error = function(e) e)
  stopifnot(
    inherits(unsupported, "error"),
    grepl("Failed to save failed fit; stopping batch", conditionMessage(unsupported))
  )
})

message("Failed-fit persistence checks passed.")
