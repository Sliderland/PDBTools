# Run from any directory with:
#   Rscript PDBTools_v3/tests/test_unconstrained_parameter_counts.R

.test_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.test_file) || !nzchar(.test_file)) {
  .test_file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(.test_file_arg)) {
    .test_file <- sub("^--file=", "", .test_file_arg[[1L]])
  }
}
.v3_dir <- dirname(normalizePath(.test_file, mustWork = TRUE))
source(file.path(.v3_dir, "..", "PDBEntryBuilder_v3.R"))

local({
  model_code <- paste(
    "parameters {",
    "  real<lower=0> sigma;",
    "  simplex[3] theta;",
    "  vector[2] beta;",
    "}",
    "transformed parameters { real beta_sum = sum(beta); }",
    "model {",
    "  sigma ~ normal(0, 1);",
    "  theta ~ dirichlet(rep_vector(1, 3));",
    "  beta ~ normal(0, 1);",
    "}",
    "generated quantities { real check_value = beta_sum; }",
    sep = "\n"
  )

  BuilderProbe <- R6::R6Class(
    "V3CountBuilderProbe",
    inherit = PDBEntryBuilder_v3,
    public = list(
      initialize = function() {
        self$stan_backend <- "rstan"
        invisible(self)
      }
    )
  )
  builder <- BuilderProbe$new()
  fit <- suppressWarnings(rstan::stan(
    model_code = model_code,
    data = list(),
    chains = 0,
    refresh = 0
  ))

  counts <- builder$generate_model_unconstrained_parameter_counts_from_fit(
    fit,
    exclude = c("beta_sum", "check_value")
  )
  stopifnot(
    identical(as.integer(counts$sigma), 1L),
    identical(as.integer(counts$theta), 2L),
    identical(as.integer(counts$beta), 2L),
    sum(unlist(counts)) == rstan::get_num_upars(fit),
    !"lp__" %in% names(counts)
  )

  selected <- builder$generate_model_unconstrained_parameter_counts_from_fit(
    fit,
    include = c("theta", "sigma")
  )
  stopifnot(
    identical(names(selected), c("theta", "sigma")),
    identical(as.integer(unlist(selected)), c(2L, 1L))
  )

  inferred <- builder$generate_model_unconstrained_parameter_counts(
    model_code = model_code,
    data = list(),
    include = c("theta", "sigma")
  )
  stopifnot(identical(as.integer(unlist(inferred)), c(2L, 1L)))

  stan_file <- tempfile(fileext = ".stan")
  writeLines(model_code, stan_file)
  stopifnot(identical(
    as.integer(unlist(builder$generate_model_unconstrained_parameter_counts(
      model_code = stan_file,
      data = list(),
      include = c("theta", "sigma")
    ))),
    c(2L, 1L)
  ))
  unlink(stan_file)

  cmdstan_available <- requireNamespace("cmdstanr", quietly = TRUE) &&
    tryCatch({
      cmdstanr::cmdstan_path()
      TRUE
    }, error = function(error) FALSE)
  if (cmdstan_available) {
    cmdstan_counts <- builder$generate_model_unconstrained_parameter_counts(
      model_code = model_code,
      data = list(),
      exclude = c("beta_sum", "check_value"),
      backend = "cmdstanr"
    )
    stopifnot(
      identical(as.integer(cmdstan_counts$sigma), 1L),
      identical(as.integer(cmdstan_counts$theta), 2L),
      identical(as.integer(cmdstan_counts$beta), 2L)
    )
  }

  expected <- list(theta = 2L, sigma = 1L)
  CountCheckProbe <- R6::R6Class(
    "V3CountCheckProbe",
    inherit = PDBEntryBuilder_v3,
    public = list(
      initialize = function() invisible(self),
      generate_model_unconstrained_parameter_counts_from_fit = function(
        fit,
        include = NULL,
        exclude = NULL
      ) list(theta = 2L, sigma = 1L)
    )
  )
  checker <- CountCheckProbe$new()
  fake_fit <- structure(list(), model_name = "toy")
  stopifnot(isTRUE(checker$check_unconstrained_parameter_counts(
    fake_fit,
    expected_dimensions = expected
  )))
  stopifnot(!checker$check_unconstrained_parameter_counts(
    fake_fit,
    expected_dimensions = list(theta = 3L, sigma = 1L)
  ))
  stopifnot(!checker$check_unconstrained_parameter_counts(
    fake_fit,
    expected_dimensions = list(theta = 2L, absent = 1L)
  ))

  message("V3 unconstrained parameter count tests passed.")
})
