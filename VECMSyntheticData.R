# Reproducible observed data for the two PosteriorDB VECM formulations.
# The first four rows are fixed conditioning values; the likelihood starts
# at row five. The process has five variables, rank two, and three common
# stochastic trends. It is simulated in transitory form; the long-run form
# has the same likelihood with gamma_long_run[j] = gamma_transitory[j] + Pi.

generate_vecm_synthetic_data <- function(
    n_obs = 160L,
    seed = 20260917L
) {
    N <- 5L
    p <- 4L
    h <- 2L
    if (
        length(n_obs) != 1L || is.na(n_obs) ||
        n_obs != as.integer(n_obs) || n_obs <= p ||
        length(seed) != 1L || is.na(seed) ||
        seed != as.integer(seed)
    ) {
        stop("`n_obs` and `seed` must be valid integers; `n_obs` must exceed 4.")
    }
    n_obs <- as.integer(n_obs)
    seed <- as.integer(seed)

    B <- matrix(
        c(0.55, -0.10, 0.25, 0.35, -0.15, 0.20),
        nrow = N - h,
        byrow = TRUE
    )
    beta <- cbind(diag(h), -t(B))
    alpha <- rbind(
        c(-0.30, 0.01),
        c(0.02, -0.27),
        c(0.03, -0.01),
        c(-0.02, 0.02),
        c(0.01, 0.02)
    )
    Pi <- alpha %*% beta
    gamma_transitory <- lapply(c(0.18, 0.07, -0.02), function(x) {
        diag(rep(x, N))
    })
    gamma_long_run <- lapply(gamma_transitory, function(x) x + Pi)

    innovation_chol <- diag(c(0.12, 0.11, 0.10, 0.09, 0.08))
    innovation_chol[2, 1] <- 0.025
    innovation_chol[3, 1] <- -0.015
    innovation_chol[4, 2] <- 0.020
    innovation_chol[5, 3] <- 0.018
    Sigma <- innovation_chol %*% t(innovation_chol)

    identity <- diag(N)
    A <- list(
        identity + Pi + gamma_transitory[[1L]],
        gamma_transitory[[2L]] - gamma_transitory[[1L]],
        gamma_transitory[[3L]] - gamma_transitory[[2L]],
        -gamma_transitory[[3L]]
    )
    companion <- matrix(0.0, N * p, N * p)
    companion[seq_len(N), ] <- do.call(cbind, A)
    companion[(N + 1L):(N * p), seq_len(N * (p - 1L))] <-
        diag(N * (p - 1L))
    roots <- eigen(companion, only.values = TRUE)$values
    unit_roots <- Mod(roots - 1) < 1e-8
    max_stable_root <- max(Mod(roots[!unit_roots]))
    if (
        qr(Pi)$rank != h || sum(unit_roots) != N - h ||
        max_stable_root >= 1 - 1e-6
    ) {
        stop("The selected VECM parameters do not have the intended rank and roots.")
    }

    prior_rng_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
        prior_seed <- get(".Random.seed", envir = .GlobalEnv)
    }
    on.exit({
        do.call(RNGkind, as.list(prior_rng_kind))
        if (had_seed) {
            assign(".Random.seed", prior_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)

    innovations <- matrix(rnorm(n_obs * N), n_obs, N) %*%
        t(innovation_chol)
    y <- matrix(
        0.0,
        nrow = n_obs,
        ncol = N,
        dimnames = list(NULL, paste0("y", seq_len(N)))
    )
    equivalence_error <- 0.0
    innovation_error <- 0.0
    for (time in (p + 1L):n_obs) {
        short_run_transitory <- rep(0.0, N)
        short_run_long_run <- rep(0.0, N)
        for (lag in seq_len(p - 1L)) {
            lagged_difference <- y[time - lag, ] - y[time - lag - 1L, ]
            short_run_transitory <- short_run_transitory +
                drop(gamma_transitory[[lag]] %*% lagged_difference)
            short_run_long_run <- short_run_long_run +
                drop(gamma_long_run[[lag]] %*% lagged_difference)
        }
        transitory_mean <- drop(Pi %*% y[time - 1L, ]) +
            short_run_transitory
        long_run_mean <- drop(Pi %*% y[time - p, ]) + short_run_long_run
        equivalence_error <- max(
            equivalence_error,
            abs(transitory_mean - long_run_mean)
        )
        y[time, ] <- y[time - 1L, ] + transitory_mean + innovations[time, ]
        innovation_error <- max(
            innovation_error,
            abs(y[time, ] - y[time - 1L, ] -
                transitory_mean - innovations[time, ])
        )
    }
    if (
        !all(is.finite(y)) || equivalence_error > 1e-10 ||
        innovation_error > 1e-10
    ) {
        stop("The simulated VECM data failed finite or equation checks.")
    }

    list(
        stan_data = list(t = n_obs, N = N, p = p, h = h, y = y),
        truth = list(
            seed = seed,
            alpha = alpha,
            beta = beta,
            Pi = Pi,
            gamma_transitory = gamma_transitory,
            gamma_long_run = gamma_long_run,
            Sigma = Sigma
        ),
        checks = list(
            rank = qr(Pi)$rank,
            unit_root_count = sum(unit_roots),
            max_other_root_modulus = max_stable_root,
            transitory_long_run_error = equivalence_error,
            innovation_error = innovation_error
        )
    )
}
