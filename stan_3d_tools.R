parameterInput <- function(inputId, label) {
    shiny::selectizeInput(
        inputId = inputId,
        label = label,
        choices = NULL,
        selected = NULL,
        multiple = FALSE,
        options = list(
            placeholder = "Type to search parameters...",
            closeAfterSelect = TRUE,
            openOnFocus = TRUE
        )
    )
}

get_draws_df <- function(fit, inc_warmup = FALSE) {
    if (inherits(fit, "CmdStanMCMC")) {
        posterior::as_draws_df(
            fit$draws(inc_warmup = inc_warmup)
        )
    } else if (inherits(fit, "stanfit")) {
        posterior::as_draws_df(fit)
    } else {
        stop("Unsupported fit object.")
    }
}

get_diagnostics_df <- function(fit, inc_warmup = FALSE) {
    if (inherits(fit, "CmdStanMCMC")) {
        as.data.frame(
            fit$sampler_diagnostics(format = "df")
        ) |>
            select(.chain, .iteration, divergent__)
    } else if (inherits(fit, "stanfit")) {
        posterior::as_draws_df(rstan::get_sampler_params(
            fit,
            inc_warmup = FALSE
        )) |>
            select(.chain, .iteration, divergent__)
    } else {
        stop("Unsupported fit object.")
    }
}


# Build a Plotly-compatible 3D kernel-density isosurface.  The density is
# deliberately computed only when a caller requests it, because KDE is much
# more expensive than drawing the posterior points.
compute_3d_density <- function(
    draw_data,
    x_variable,
    y_variable,
    z_variable,
    grid_size = 25L,
    max_kde_points = 20000L
) {
    if (!requireNamespace("ks", quietly = TRUE)) {
        stop("Package `ks` is required for 3D density plots.")
    }

    coordinates <- draw_data[, c(x_variable, y_variable, z_variable),
                             drop = FALSE]
    coordinates <- coordinates[
        apply(coordinates, 1L, function(row) all(is.finite(row))),
        ,
        drop = FALSE
    ]

    if (nrow(coordinates) < 20L) {
        stop("At least 20 finite draws are required for a 3D density plot.")
    }

    if (any(vapply(coordinates, function(x) diff(range(x)) == 0, logical(1)))) {
        stop("A 3D density cannot be computed when an axis is constant.")
    }

    if (nrow(coordinates) > max_kde_points) {
        coordinates <- coordinates[
            sample.int(nrow(coordinates), max_kde_points),
            ,
            drop = FALSE
        ]
    }

    grid_size <- as.integer(grid_size)
    density <- ks::kde(
        x = as.matrix(coordinates),
        gridsize = rep(grid_size, 3L),
        binned = TRUE
    )

    values <- as.vector(density$estimate)
    grid <- expand.grid(density$eval.points)
    finite <- is.finite(values) & values > 0

    if (sum(finite) < 10L) {
        stop("The 3D density estimate did not contain enough finite values.")
    }

    values <- values[finite]
    grid <- grid[finite, , drop = FALSE]
    list(
        x = grid[[1L]],
        y = grid[[2L]],
        z = grid[[3L]],
        value = values,
        isomin = as.numeric(stats::quantile(values, 0.65)),
        isomax = max(values),
        surface_count = 4L
    )
}


add_density_trace <- function(plot, density, name, color) {
    plot |>
        plotly::add_trace(
            x = density$x,
            y = density$y,
            z = density$z,
            value = density$value,
            type = "isosurface",
            name = paste(name, "density"),
            isomin = density$isomin,
            isomax = density$isomax,
            surface = list(count = density$surface_count),
            opacity = 0.35,
            colorscale = list(
                list(0, color),
                list(1, color)
            ),
            showscale = FALSE,
            hoverinfo = "skip"
        )
}


get_cached_3d_density <- function(
    cache,
    cache_key,
    draw_data,
    x_variable,
    y_variable,
    z_variable,
    grid_size
) {
    if (exists(cache_key, envir = cache, inherits = FALSE)) {
        return(get(cache_key, envir = cache, inherits = FALSE))
    }

    density <- tryCatch(
        compute_3d_density(
            draw_data,
            x_variable,
            y_variable,
            z_variable,
            grid_size = grid_size
        ),
        error = function(error) error
    )
    assign(cache_key, density, envir = cache)
    density
}


# Read a local PosteriorDB reference-draw archive without requiring its
# diagnostic metadata.  This is useful for older/local archives whose
# metadata predates fields expected by the current posteriordb package.
read_pdb_reference_draws <- function(
    pdb_path,
    reference_posterior_name
) {
    pdb_path <- path.expand(pdb_path)
    if (!dir.exists(pdb_path)) {
        stop("PosteriorDB path does not exist: ", pdb_path)
    }

    if (
        length(reference_posterior_name) != 1L ||
            is.na(reference_posterior_name) ||
            !nzchar(reference_posterior_name)
    ) {
        stop("`reference_posterior_name` must be one non-empty string.")
    }

    archive_path <- file.path(
        pdb_path,
        "posterior_database",
        "reference_posteriors",
        "draws",
        "draws",
        paste0(reference_posterior_name, ".json.zip")
    )
    info_path <- file.path(
        pdb_path,
        "posterior_database",
        "reference_posteriors",
        "draws",
        "info",
        paste0(reference_posterior_name, ".info.json")
    )

    if (!file.exists(archive_path)) {
        stop(
            "Reference-draw archive not found for `",
            reference_posterior_name,
            "`: ",
            archive_path
        )
    }

    if (!file.exists(info_path)) {
        stop(
            "Reference-draw metadata not found for `",
            reference_posterior_name,
            "`: ",
            info_path
        )
    }

    reference_info <- jsonlite::read_json(
        info_path,
        simplifyVector = FALSE
    )
    variable_names <- reference_info$diagnostics$
        diagnostic_information$names

    if (is.null(variable_names) || length(variable_names) == 0L) {
        stop(
            "Reference-draw metadata does not contain variable names: ",
            info_path
        )
    }

    archive_members <- unzip(archive_path, list = TRUE)$Name
    json_member <- archive_members[
        grepl("\\.json$", archive_members, ignore.case = TRUE)
    ]

    if (length(json_member) != 1L) {
        stop(
            "Expected exactly one JSON member in reference-draw archive: ",
            archive_path
        )
    }

    extraction_dir <- tempfile("pdb_reference_draws_")
    dir.create(extraction_dir)
    on.exit(unlink(extraction_dir, recursive = TRUE), add = TRUE)

    unzip(
        archive_path,
        files = json_member,
        exdir = extraction_dir,
        junkpaths = TRUE
    )
    json_path <- file.path(extraction_dir, basename(json_member))

    raw_draws <- jsonlite::read_json(
        json_path,
        simplifyVector = FALSE
    )

    dimensions <- c(
        length(raw_draws),
        length(raw_draws[[1L]]),
        length(raw_draws[[1L]][[1L]])
    )
    number_of_chains <- reference_info$inference$
        method_arguments$chains
    number_of_chains <- as.integer(number_of_chains)

    if (
        length(number_of_chains) != 1L ||
            is.na(number_of_chains) ||
            number_of_chains < 1L
    ) {
        stop("Reference-draw metadata does not contain a valid chain count.")
    }

    if (
        dimensions[1L] == number_of_chains &&
            dimensions[2L] == length(variable_names)
    ) {
        # Standard PosteriorDB layout: chain x variable x iteration.
        values <- aperm(
            array(
                unlist(raw_draws, use.names = FALSE),
                dim = c(
                    dimensions[3L],
                    dimensions[2L],
                    dimensions[1L]
                )
            ),
            perm = c(1L, 3L, 2L)
        )
    } else if (
        dimensions[2L] == number_of_chains &&
            dimensions[3L] == length(variable_names)
    ) {
        # Heaps archives in this local database: iteration x chain x variable.
        values <- aperm(
            array(
                unlist(raw_draws, use.names = FALSE),
                dim = dimensions[c(3L, 2L, 1L)]
            ),
            perm = c(3L, 2L, 1L)
        )
    } else {
        stop(
            "Could not reconcile draw archive dimensions ",
            paste(dimensions, collapse = " x "),
            " with ",
            number_of_chains,
            " chains and ",
            length(variable_names),
            " variables."
        )
    }

    dimnames(values) <- list(
        iteration = as.character(seq_len(dim(values)[1L])),
        chain = as.character(seq_len(dim(values)[2L])),
        variable = variable_names
    )

    posterior::as_draws_array(values)
}


launch_pdb_reference_3d <- function(
    pdb_path,
    reference_posterior_name,
    max_points = 5000
) {
    draws <- read_pdb_reference_draws(
        pdb_path = pdb_path,
        reference_posterior_name = reference_posterior_name
    )

    launch_draws_3d(
        draws = draws,
        max_points = max_points
    )
}


read_pdb_github_reference_draws <- function(
    pdb,
    reference_posterior_name
) {
    if (!inherits(pdb, "pdb_github")) {
        stop("`pdb` must be a posteriordb `pdb_github` connection.")
    }

    temporary_root <- tempfile("pdb_github_reference_")
    draw_dir <- file.path(
        temporary_root,
        "posterior_database",
        "reference_posteriors",
        "draws",
        "draws"
    )
    info_dir <- file.path(
        temporary_root,
        "posterior_database",
        "reference_posteriors",
        "draws",
        "info"
    )
    dir.create(draw_dir, recursive = TRUE)
    dir.create(info_dir, recursive = TRUE)
    on.exit(unlink(temporary_root, recursive = TRUE), add = TRUE)

    archive_name <- paste0(reference_posterior_name, ".json.zip")
    info_name <- paste0(reference_posterior_name, ".info.json")
    archive_path <- file.path(draw_dir, archive_name)
    info_path <- file.path(info_dir, info_name)

    # Use posteriordb's persistent connection cache so repeated launches do
    # not repeatedly call the GitHub API or redownload the same large archive.
    archive_cached_path <- posteriordb:::pdb_cached_local_file_path(
        pdb,
        file.path(
            "reference_posteriors",
            "draws",
            "draws",
            archive_name
        ),
        unzip = FALSE
    )
    info_cached_path <- posteriordb:::pdb_cached_local_file_path(
        pdb,
        file.path(
            "reference_posteriors",
            "draws",
            "info",
            info_name
        ),
        unzip = FALSE
    )
    archive_downloaded <- file.copy(
        archive_cached_path,
        archive_path,
        overwrite = TRUE
    )
    info_downloaded <- file.copy(
        info_cached_path,
        info_path,
        overwrite = TRUE
    )

    if (!isTRUE(archive_downloaded) || !file.exists(archive_path)) {
        stop(
            "Could not download reference-draw archive for `",
            reference_posterior_name,
            "` from the GitHub PosteriorDB connection."
        )
    }
    if (!isTRUE(info_downloaded) || !file.exists(info_path)) {
        stop(
            "Could not download reference-draw metadata for `",
            reference_posterior_name,
            "` from the GitHub PosteriorDB connection."
        )
    }

    read_pdb_reference_draws(
        pdb_path = temporary_root,
        reference_posterior_name = reference_posterior_name
    )
}


launch_pdb_github_reference_3d <- function(
    pdb,
    reference_posterior_name,
    max_points = 5000
) {
    draws <- read_pdb_github_reference_draws(
        pdb = pdb,
        reference_posterior_name = reference_posterior_name
    )

    launch_draws_3d(
        draws = draws,
        max_points = max_points
    )
}


launch_stan_3d <- function(fit) {
    library(shiny)
    library(plotly)
    library(dplyr)

    # Convert posterior's draws_df to an ordinary data frame before joining.
    # This avoids warnings caused by removing draws_df metadata.
    draws <- get_draws_df(fit)
    diagnostics <- get_diagnostics_df(fit)

    plot_data <- draws |>
        left_join(
            diagnostics,
            by = c(".chain", ".iteration")
        ) |>
        mutate(
            divergence = factor(
                divergent__,
                levels = c(0, 1),
                labels = c("Regular", "Divergent")
            )
        )

    excluded <- c(
        ".chain",
        ".iteration",
        ".draw",
        "lp__",
        "divergent__",
        "divergence"
    )

    parameters <- names(plot_data)[
        !names(plot_data) %in% excluded &
            !grepl("__$", names(plot_data))
    ]

    if (length(parameters) < 3) {
        stop("The fit must contain at least three selectable variables.")
    }

    ui <- fluidPage(
        titlePanel("Interactive Stan posterior explorer"),

        sidebarLayout(
            sidebarPanel(
                parameterInput(
                    "x",
                    "X parameter"
                ),
                parameterInput(
                    "y",
                    "Y parameter"
                ),
                parameterInput(
                    "z",
                    "Z parameter"
                ),

                radioButtons(
                    "colour",
                    "Colour points by",
                    choices = c("Divergence", "Chain"),
                    selected = if (
                        any(plot_data$divergent__ == 1, na.rm = TRUE)
                    ) {
                        "Divergence"
                    } else {
                        "Chain"
                    }
                ),

                radioButtons(
                    "display_mode",
                    "Display",
                    choices = c("Scatter", "Density", "Both"),
                    selected = "Scatter"
                ),

                sliderInput(
                    "density_grid",
                    "Density grid resolution",
                    min = 15,
                    max = 40,
                    value = 25,
                    step = 1
                ),

                sliderInput(
                    "opacity",
                    "Point opacity",
                    min = 0.05,
                    max = 1,
                    value = 0.6
                ),

                sliderInput(
                    "size",
                    "Point size",
                    min = 1,
                    max = 8,
                    value = 2
                )
            ),

            mainPanel(
                plotlyOutput("posterior_plot", height = "750px")
            )
        )
    )

    server <- function(input, output, session) {
        session$onFlushed(
            function() {
                updateSelectizeInput(
                    session,
                    "x",
                    choices = parameters,
                    selected = parameters[1],
                    server = TRUE
                )
                updateSelectizeInput(
                    session,
                    "y",
                    choices = parameters,
                    selected = parameters[2],
                    server = TRUE
                )
                updateSelectizeInput(
                    session,
                    "z",
                    choices = parameters,
                    selected = parameters[3],
                    server = TRUE
                )
            },
            once = TRUE
        )

        density_cache <- new.env(parent = emptyenv())

        output$posterior_plot <- renderPlotly({
            req(input$x, input$y, input$z)

            show_scatter <- input$display_mode %in% c("Scatter", "Both")
            show_density <- input$display_mode %in% c("Density", "Both")
            plot <- plotly::plot_ly()

            if (show_scatter) {
                colour_variable <- if (input$colour == "Divergence") {
                    plot_data$divergence
                } else {
                    factor(plot_data$.chain)
                }

                plot <- plot |>
                    plotly::add_trace(
                        x = plot_data[[input$x]],
                        y = plot_data[[input$y]],
                        z = plot_data[[input$z]],
                        color = colour_variable,
                        type = "scatter3d",
                        mode = "markers",
                        marker = list(
                            size = input$size,
                            opacity = input$opacity
                        )
                    )
            }

            if (show_density) {
                density <- get_cached_3d_density(
                    density_cache,
                    paste(input$x, input$y, input$z, input$density_grid, sep = "|"),
                    plot_data,
                    input$x,
                    input$y,
                    input$z,
                    input$density_grid
                )
                validate(
                    shiny::need(
                        !inherits(density, "error"),
                        paste("Density unavailable:", density$message)
                    )
                )
                plot <- add_density_trace(
                    plot,
                    density,
                    "Posterior",
                    "#2563EB"
                )
            }

            plot |>
                plotly::layout(
                    scene = list(
                        xaxis = list(title = input$x),
                        yaxis = list(title = input$y),
                        zaxis = list(title = input$z)
                    ),
                    legend = list(title = list(text = input$colour))
                )
        })
    }

    shinyApp(ui, server)
}


launch_fit_comparison <- function(fit1, fit2) {
    library(shiny)
    library(plotly)

    draws1 <- as.data.frame(fit1$draws(format = "df"))
    draws2 <- as.data.frame(fit2$draws(format = "df"))

    excluded <- c(".chain", ".iteration", ".draw", "lp__")

    variables1 <- names(draws1)[
        !names(draws1) %in% excluded &
            !grepl("__$", names(draws1))
    ]

    variables2 <- names(draws2)[
        !names(draws2) %in% excluded &
            !grepl("__$", names(draws2))
    ]

    if (length(variables1) < 3 || length(variables2) < 3) {
        stop("Each fit must contain at least three selectable variables.")
    }

    ui <- fluidPage(
        titlePanel("Compare two Stan posterior distributions"),

        sidebarLayout(
            sidebarPanel(
                h4("X axis"),
                fluidRow(
                    column(
                        6,
                        parameterInput(
                            "x1",
                            "Fit 1"
                        )
                    ),
                    column(
                        6,
                        parameterInput(
                            "x2",
                            "Fit 2"
                        )
                    )
                ),

                h4("Y axis"),
                fluidRow(
                    column(
                        6,
                        parameterInput(
                            "y1",
                            "Fit 1"
                        )
                    ),
                    column(
                        6,
                        parameterInput(
                            "y2",
                            "Fit 2"
                        )
                    )
                ),

                h4("Z axis"),
                fluidRow(
                    column(
                        6,
                        parameterInput(
                            "z1",
                            "Fit 1"
                        )
                    ),
                    column(
                        6,
                        parameterInput(
                            "z2",
                            "Fit 2"
                        )
                    )
                ),

                sliderInput(
                    "point_size",
                    "Point size",
                    min = 1,
                    max = 8,
                    value = 2
                ),

                sliderInput(
                    "opacity",
                    "Opacity",
                    min = 0.05,
                    max = 1,
                    value = 0.45
                ),

                checkboxInput("show_fit1", "Show Fit 1", value = TRUE),
                checkboxInput("show_fit2", "Show Fit 2", value = TRUE),
                checkboxInput(
                    "show_density1",
                    "Show Fit 1 density",
                    value = FALSE
                ),
                checkboxInput(
                    "show_density2",
                    "Show Fit 2 density",
                    value = FALSE
                ),
                sliderInput(
                    "density_grid",
                    "Density grid resolution",
                    min = 15,
                    max = 40,
                    value = 25,
                    step = 1
                )
            ),

            mainPanel(
                plotlyOutput("posterior_plot", height = "800px")
            )
        )
    )

    server <- function(input, output, session) {
        session$onFlushed(
            function() {
                initial_values1 <- variables1[1:3]
                initial_values2 <- variables2[1:3]

                for (axis_index in seq_along(c("x", "y", "z"))) {
                    axis <- c("x", "y", "z")[axis_index]

                    updateSelectizeInput(
                        session,
                        paste0(axis, "1"),
                        choices = variables1,
                        selected = initial_values1[axis_index],
                        server = TRUE
                    )
                    updateSelectizeInput(
                        session,
                        paste0(axis, "2"),
                        choices = variables2,
                        selected = initial_values2[axis_index],
                        server = TRUE
                    )
                }
            },
            once = TRUE
        )

        observeEvent(input$x1, {
            if (!is.null(input$x1) && input$x1 %in% variables2) {
                updateSelectizeInput(
                    session,
                    inputId = "x2",
                    choices = variables2,
                    selected = input$x1,
                    server = TRUE
                )
            }
        })

        observeEvent(input$y1, {
            if (!is.null(input$y1) && input$y1 %in% variables2) {
                updateSelectizeInput(
                    session,
                    inputId = "y2",
                    choices = variables2,
                    selected = input$y1,
                    server = TRUE
                )
            }
        })

        observeEvent(input$z1, {
            if (!is.null(input$z1) && input$z1 %in% variables2) {
                updateSelectizeInput(
                    session,
                    inputId = "z2",
                    choices = variables2,
                    selected = input$z1,
                    server = TRUE
                )
            }
        })

        density_cache <- new.env(parent = emptyenv())

        output$posterior_plot <- renderPlotly({
            req(input$x1, input$x2, input$y1, input$y2, input$z1, input$z2)

            plot <- plot_ly()

            if (input$show_fit1) {
                plot <- plot |>
                    add_trace(
                        x = draws1[[input$x1]],
                        y = draws1[[input$y1]],
                        z = draws1[[input$z1]],
                        type = "scatter3d",
                        mode = "markers",
                        name = "Fit 1",
                        marker = list(
                            size = input$point_size,
                            opacity = input$opacity,
                            color = "#2563EB"
                        ),
                        text = paste0(
                            input$x1,
                            ": ",
                            signif(draws1[[input$x1]], 4),
                            "<br>",
                            input$y1,
                            ": ",
                            signif(draws1[[input$y1]], 4),
                            "<br>",
                            input$z1,
                            ": ",
                            signif(draws1[[input$z1]], 4),
                            "<br>Chain: ",
                            draws1$.chain
                        ),
                        hoverinfo = "text"
                    )
            }

            if (isTRUE(input$show_density1)) {
                density1 <- get_cached_3d_density(
                    density_cache,
                    paste("fit1", input$x1, input$y1, input$z1,
                          input$density_grid, sep = "|"),
                    draws1,
                    input$x1,
                    input$y1,
                    input$z1,
                    input$density_grid
                )
                validate(
                    need(
                        !inherits(density1, "error"),
                        paste("Fit 1 density unavailable:", density1$message)
                    )
                )
                plot <- add_density_trace(
                    plot,
                    density1,
                    "Fit 1",
                    "#2563EB"
                )
            }

            if (input$show_fit2) {
                plot <- plot |>
                    add_trace(
                        x = draws2[[input$x2]],
                        y = draws2[[input$y2]],
                        z = draws2[[input$z2]],
                        type = "scatter3d",
                        mode = "markers",
                        name = "Fit 2",
                        marker = list(
                            size = input$point_size,
                            opacity = input$opacity,
                            color = "#DC2626"
                        ),
                        text = paste0(
                            input$x2,
                            ": ",
                            signif(draws2[[input$x2]], 4),
                            "<br>",
                            input$y2,
                            ": ",
                            signif(draws2[[input$y2]], 4),
                            "<br>",
                            input$z2,
                            ": ",
                            signif(draws2[[input$z2]], 4),
                            "<br>Chain: ",
                            draws2$.chain
                        ),
                        hoverinfo = "text"
                    )
            }

            if (isTRUE(input$show_density2)) {
                density2 <- get_cached_3d_density(
                    density_cache,
                    paste("fit2", input$x2, input$y2, input$z2,
                          input$density_grid, sep = "|"),
                    draws2,
                    input$x2,
                    input$y2,
                    input$z2,
                    input$density_grid
                )
                validate(
                    need(
                        !inherits(density2, "error"),
                        paste("Fit 2 density unavailable:", density2$message)
                    )
                )
                plot <- add_density_trace(
                    plot,
                    density2,
                    "Fit 2",
                    "#DC2626"
                )
            }

            plot |>
                layout(
                    scene = list(
                        xaxis = list(title = paste(input$x1, "↔", input$x2)),
                        yaxis = list(title = paste(input$y1, "↔", input$y2)),
                        zaxis = list(title = paste(input$z1, "↔", input$z2))
                    ),
                    legend = list(title = list(text = "Posterior"))
                )
        })
    }

    shinyApp(ui, server)
}


launch_multi_fit_comparison <- function(
    fits,
    max_points_per_fit = 5000
) {
    library(shiny)
    library(plotly)

    if (!is.list(fits) || length(fits) < 1) {
        stop("`fits` must be a non-empty list of CmdStanR fit objects.")
    }

    if (
        length(max_points_per_fit) != 1 ||
            is.na(max_points_per_fit) ||
            max_points_per_fit <= 0
    ) {
        stop("`max_points_per_fit` must be a positive number or Inf.")
    }

    number_of_fits <- length(fits)
    fit_names <- names(fits)

    if (is.null(fit_names)) {
        fit_names <- rep("", number_of_fits)
    }

    missing_names <- is.na(fit_names) | fit_names == ""
    fit_names[missing_names] <- paste("Fit", which(missing_names))

    draws <- lapply(seq_along(fits), function(i) {
        tryCatch(
            as.data.frame(fits[[i]]$draws(format = "df")),
            error = function(error) {
                stop(
                    "Could not extract draws from ",
                    fit_names[i],
                    ": ",
                    conditionMessage(error)
                )
            }
        )
    })

    excluded <- c(".chain", ".iteration", ".draw", "lp__")

    variables <- lapply(draws, function(draw_data) {
        column_names <- names(draw_data)

        column_names[
            !column_names %in% excluded &
                !grepl("__$", column_names)
        ]
    })

    too_small <- which(lengths(variables) < 3)

    if (length(too_small) > 0) {
        stop(
            "Each fit must contain at least three selectable variables. ",
            "Too few were found in: ",
            paste(fit_names[too_small], collapse = ", ")
        )
    }

    # Downsample only the points sent to Plotly. The original fit objects and
    # extracted draws remain unchanged.
    plot_draws <- lapply(draws, function(draw_data) {
        if (
            is.finite(max_points_per_fit) &&
                nrow(draw_data) > max_points_per_fit
        ) {
            indices <- unique(
                round(
                    seq(
                        1,
                        nrow(draw_data),
                        length.out = max_points_per_fit
                    )
                )
            )

            draw_data[indices, , drop = FALSE]
        } else {
            draw_data
        }
    })

    fit_colors <- grDevices::hcl.colors(
        number_of_fits,
        palette = "Dark 3"
    )

    axes <- c("x", "y", "z")

    ui <- fluidPage(
        titlePanel("Compare multiple Stan posterior distributions"),

        sidebarLayout(
            sidebarPanel(
                width = 4,

                helpText(
                    "The first fit is the reference. Exact parameter-name",
                    "matches are selected automatically in the other fits."
                ),

                lapply(seq_len(number_of_fits), function(i) {
                    wellPanel(
                        h4(
                            fit_names[i],
                            style = paste0("color: ", fit_colors[i], ";")
                        ),

                        fluidRow(
                            column(
                                4,
                                parameterInput(
                                    paste0("x_", i),
                                    "X"
                                )
                            ),
                            column(
                                4,
                                parameterInput(
                                    paste0("y_", i),
                                    "Y"
                                )
                            ),
                            column(
                                4,
                                parameterInput(
                                    paste0("z_", i),
                                    "Z"
                                )
                            )
                        ),

                        checkboxInput(
                            paste0("show_", i),
                            paste("Show", fit_names[i]),
                            value = TRUE
                        ),

                        checkboxInput(
                            paste0("show_density_", i),
                            paste("Show density", fit_names[i]),
                            value = FALSE
                        )
                    )
                }),

                sliderInput(
                    "point_size",
                    "Point size",
                    min = 1,
                    max = 8,
                    value = 2
                ),

                sliderInput(
                    "opacity",
                    "Opacity",
                    min = 0.05,
                    max = 1,
                    value = 0.45
                ),

                sliderInput(
                    "density_grid",
                    "Density grid resolution",
                    min = 15,
                    max = 40,
                    value = 25,
                    step = 1
                )
            ),

            mainPanel(
                width = 8,
                plotlyOutput("posterior_plot", height = "850px")
            )
        )
    )

    server <- function(input, output, session) {
        session$onFlushed(
            function() {
                reference_initial <- variables[[1]][1:3]

                for (i in seq_len(number_of_fits)) {
                    for (axis_index in seq_along(axes)) {
                        reference_value <- reference_initial[axis_index]

                        selected_value <- if (
                            reference_value %in% variables[[i]]
                        ) {
                            reference_value
                        } else {
                            variables[[i]][axis_index]
                        }

                        updateSelectizeInput(
                            session,
                            inputId = paste0(axes[axis_index], "_", i),
                            choices = variables[[i]],
                            selected = selected_value,
                            server = TRUE
                        )
                    }
                }
            },
            once = TRUE
        )

        # Changes to the reference fit propagate only when an exact match is
        # available. Otherwise, the other fit's selection is left untouched.
        lapply(axes, function(axis) {
            observeEvent(
                input[[paste0(axis, "_1")]],
                {
                    selected_value <- input[[paste0(axis, "_1")]]

                    if (!is.null(selected_value) && number_of_fits > 1) {
                        for (i in 2:number_of_fits) {
                            if (selected_value %in% variables[[i]]) {
                                updateSelectizeInput(
                                    session,
                                    inputId = paste0(axis, "_", i),
                                    choices = variables[[i]],
                                    selected = selected_value,
                                    server = TRUE
                                )
                            }
                        }
                    }
                },
                ignoreInit = TRUE
            )
        })

        density_cache <- new.env(parent = emptyenv())

        output$posterior_plot <- renderPlotly({
            plot <- plot_ly()
            visible_fits <- 0L

            for (i in seq_len(number_of_fits)) {
                show_scatter <- isTRUE(input[[paste0("show_", i)]])
                show_density <- isTRUE(
                    input[[paste0("show_density_", i)]]
                )

                if (!show_scatter && !show_density) {
                    next
                }

                x_variable <- input[[paste0("x_", i)]]
                y_variable <- input[[paste0("y_", i)]]
                z_variable <- input[[paste0("z_", i)]]

                req(x_variable, y_variable, z_variable)

                draw_data <- plot_draws[[i]]

                hover_text <- paste0(
                    "<b>",
                    fit_names[i],
                    "</b>",
                    "<br>",
                    x_variable,
                    ": ",
                    signif(draw_data[[x_variable]], 5),
                    "<br>",
                    y_variable,
                    ": ",
                    signif(draw_data[[y_variable]], 5),
                    "<br>",
                    z_variable,
                    ": ",
                    signif(draw_data[[z_variable]], 5),
                    "<br>Chain: ",
                    draw_data$.chain,
                    "<br>Iteration: ",
                    draw_data$.iteration
                )

                if (show_scatter) {
                    plot <- plot |>
                        add_trace(
                            x = draw_data[[x_variable]],
                            y = draw_data[[y_variable]],
                            z = draw_data[[z_variable]],
                            type = "scatter3d",
                            mode = "markers",
                            name = fit_names[i],
                            marker = list(
                                size = input$point_size,
                                opacity = input$opacity,
                                color = fit_colors[i]
                            ),
                            text = hover_text,
                            hoverinfo = "text"
                        )
                }

                if (show_density) {
                    density <- get_cached_3d_density(
                        density_cache,
                        paste(
                            i,
                            x_variable,
                            y_variable,
                            z_variable,
                            input$density_grid,
                            sep = "|"
                        ),
                        draw_data,
                        x_variable,
                        y_variable,
                        z_variable,
                        input$density_grid
                    )
                    validate(
                        need(
                            !inherits(density, "error"),
                            paste(
                                fit_names[i],
                                "density unavailable:",
                                density$message
                            )
                        )
                    )
                    plot <- add_density_trace(
                        plot,
                        density,
                        fit_names[i],
                        fit_colors[i]
                    )
                }

                visible_fits <- visible_fits + 1L
            }

            validate(
                need(visible_fits > 0, "Select at least one fit to display.")
            )

            plot |>
                layout(
                    scene = list(
                        xaxis = list(title = "X"),
                        yaxis = list(title = "Y"),
                        zaxis = list(title = "Z")
                    ),
                    legend = list(title = list(text = "Posterior"))
                )
        })
    }

    shinyApp(ui, server)
}


# Explore any posterior::draws object without first materializing every
# parameter as a data frame. This is particularly useful for wide PosteriorDB
# reference-draw objects.
launch_draws_3d <- function(
    draws,
    max_points = 5000
) {
    if (!posterior::is_draws(draws)) {
        stop("`draws` must be a posterior::draws object.")
    }

    if (
        length(max_points) != 1 ||
            is.na(max_points) ||
            max_points <= 0
    ) {
        stop("`max_points` must be a positive number or Inf.")
    }

    draw_variables <- posterior::variables(draws)
    excluded <- c(
        ".chain",
        ".iteration",
        ".draw",
        "lp__",
        "divergent__"
    )

    parameters <- draw_variables[
        !draw_variables %in% excluded &
            !grepl("__$", draw_variables)
    ]

    if (length(parameters) < 3) {
        stop("`draws` must contain at least three selectable variables.")
    }

    has_divergences <- "divergent__" %in% draw_variables
    colour_choices <- "Chain"

    if (has_divergences) {
        colour_choices <- c("Divergence", "Chain")
    }

    ui <- shiny::fluidPage(
        shiny::titlePanel("Interactive posterior draws explorer"),

        shiny::sidebarLayout(
            shiny::sidebarPanel(
                parameterInput("draws_x", "X parameter"),
                parameterInput("draws_y", "Y parameter"),
                parameterInput("draws_z", "Z parameter"),

                shiny::radioButtons(
                    "draws_colour",
                    "Colour points by",
                    choices = colour_choices,
                    selected = colour_choices[1]
                ),

                shiny::radioButtons(
                    "draws_display_mode",
                    "Display",
                    choices = c("Scatter", "Density", "Both"),
                    selected = "Scatter"
                ),

                shiny::sliderInput(
                    "draws_density_grid",
                    "Density grid resolution",
                    min = 15,
                    max = 40,
                    value = 25,
                    step = 1
                ),

                shiny::sliderInput(
                    "draws_opacity",
                    "Point opacity",
                    min = 0.05,
                    max = 1,
                    value = 0.6
                ),

                shiny::sliderInput(
                    "draws_size",
                    "Point size",
                    min = 1,
                    max = 8,
                    value = 2
                ),

                shiny::helpText(
                    paste0(
                        posterior::ndraws(draws),
                        " draws in ",
                        posterior::nchains(draws),
                        " chains; at most ",
                        if (is.finite(max_points)) max_points else "all",
                        " points are plotted."
                    )
                ),

                if (!has_divergences) {
                    shiny::helpText(
                        paste(
                            "Per-draw divergent__ diagnostics are not present;",
                            "divergence colouring is unavailable."
                        )
                    )
                }
            ),

            shiny::mainPanel(
                plotly::plotlyOutput("draws_posterior_plot", height = "750px")
            )
        )
    )

    server <- function(input, output, session) {
        session$onFlushed(
            function() {
                shiny::updateSelectizeInput(
                    session,
                    "draws_x",
                    choices = parameters,
                    selected = parameters[1],
                    server = TRUE
                )
                shiny::updateSelectizeInput(
                    session,
                    "draws_y",
                    choices = parameters,
                    selected = parameters[2],
                    server = TRUE
                )
                shiny::updateSelectizeInput(
                    session,
                    "draws_z",
                    choices = parameters,
                    selected = parameters[3],
                    server = TRUE
                )
            },
            once = TRUE
        )

        selected_draws <- shiny::reactive({
            shiny::req(input$draws_x, input$draws_y, input$draws_z)

            selected_variables <- unique(
                c(
                    input$draws_x,
                    input$draws_y,
                    input$draws_z,
                    if (has_divergences) "divergent__"
                )
            )

            plot_data <- draws |>
                posterior::subset_draws(variable = selected_variables) |>
                posterior::as_draws_df() |>
                as.data.frame()

            if (
                is.finite(max_points) &&
                    nrow(plot_data) > max_points
            ) {
                # When available, retain divergent transitions preferentially
                # so rare divergences are not hidden by display downsampling.
                divergent_indices <- integer()

                if (has_divergences) {
                    divergent_indices <- which(plot_data$divergent__ == 1)
                }

                if (length(divergent_indices) >= max_points) {
                    keep <- divergent_indices[
                        unique(
                            round(
                                seq(
                                    1,
                                    length(divergent_indices),
                                    length.out = max_points
                                )
                            )
                        )
                    ]
                } else {
                    regular_indices <- setdiff(
                        seq_len(nrow(plot_data)),
                        divergent_indices
                    )
                    number_regular <- max_points - length(divergent_indices)
                    regular_keep <- regular_indices[
                        unique(
                            round(
                                seq(
                                    1,
                                    length(regular_indices),
                                    length.out = number_regular
                                )
                            )
                        )
                    ]
                    keep <- sort(c(divergent_indices, regular_keep))
                }

                plot_data <- plot_data[keep, , drop = FALSE]
            }

            plot_data
        })

        density_cache <- new.env(parent = emptyenv())

        output$draws_posterior_plot <- plotly::renderPlotly({
            plot_data <- selected_draws()

            show_scatter <- input$draws_display_mode %in% c(
                "Scatter",
                "Both"
            )
            show_density <- input$draws_display_mode %in% c(
                "Density",
                "Both"
            )
            plot <- plotly::plot_ly()

            if (show_scatter) {
                colour_variable <- if (
                    has_divergences &&
                        identical(input$draws_colour, "Divergence")
                ) {
                    factor(
                        plot_data$divergent__,
                        levels = c(0, 1),
                        labels = c("Regular", "Divergent")
                    )
                } else {
                    factor(plot_data$.chain)
                }

                hover_text <- paste0(
                    input$draws_x,
                    ": ",
                    signif(plot_data[[input$draws_x]], 4),
                    "<br>",
                    input$draws_y,
                    ": ",
                    signif(plot_data[[input$draws_y]], 4),
                    "<br>",
                    input$draws_z,
                    ": ",
                    signif(plot_data[[input$draws_z]], 4),
                    "<br>Chain: ",
                    plot_data$.chain,
                    "<br>Iteration: ",
                    plot_data$.iteration
                )

                plot <- plot |>
                    plotly::add_trace(
                        x = plot_data[[input$draws_x]],
                        y = plot_data[[input$draws_y]],
                        z = plot_data[[input$draws_z]],
                        color = colour_variable,
                        type = "scatter3d",
                        mode = "markers",
                        marker = list(
                            size = input$draws_size,
                            opacity = input$draws_opacity
                        ),
                        text = hover_text,
                        hoverinfo = "text"
                    )
            }

            if (show_density) {
                density <- get_cached_3d_density(
                    density_cache,
                    paste(
                        input$draws_x,
                        input$draws_y,
                        input$draws_z,
                        input$draws_density_grid,
                        sep = "|"
                    ),
                    plot_data,
                    input$draws_x,
                    input$draws_y,
                    input$draws_z,
                    input$draws_density_grid
                )
                shiny::validate(
                    shiny::need(
                        !inherits(density, "error"),
                        paste("Density unavailable:", density$message)
                    )
                )
                plot <- add_density_trace(
                    plot,
                    density,
                    "Posterior",
                    "#2563EB"
                )
            }

            plot |>
                plotly::layout(
                    scene = list(
                        xaxis = list(title = input$draws_x),
                        yaxis = list(title = input$draws_y),
                        zaxis = list(title = input$draws_z)
                    ),
                    legend = list(
                        title = list(text = input$draws_colour)
                    )
                )
        })
    }

    shiny::shinyApp(ui, server)
}

# parameterInput <- function(inputId, label) {
#     shiny::selectizeInput(
#         inputId = inputId,
#         label = label,
#         choices = NULL,
#         selected = NULL,
#         multiple = FALSE,
#         options = list(
#             placeholder = "Type to search parameters...",
#             closeAfterSelect = TRUE,
#             openOnFocus = TRUE
#         )
#     )
# }

# launch_stan_3d <- function(fit) {
#     library(shiny)
#     library(plotly)
#     library(dplyr)

#     # Convert posterior's draws_df to an ordinary data frame before joining.
#     # This avoids warnings caused by removing draws_df metadata.
#     draws <- as.data.frame(fit$draws(format = "df"))

#     diagnostics <- as.data.frame(
#         fit$sampler_diagnostics(format = "df")
#     ) |>
#         select(.chain, .iteration, divergent__)

#     plot_data <- draws |>
#         left_join(
#             diagnostics,
#             by = c(".chain", ".iteration")
#         ) |>
#         mutate(
#             divergence = factor(
#                 divergent__,
#                 levels = c(0, 1),
#                 labels = c("Regular", "Divergent")
#             )
#         )

#     excluded <- c(
#         ".chain",
#         ".iteration",
#         ".draw",
#         "lp__",
#         "divergent__",
#         "divergence"
#     )

#     parameters <- names(plot_data)[
#         !names(plot_data) %in% excluded &
#             !grepl("__$", names(plot_data))
#     ]

#     if (length(parameters) < 3) {
#         stop("The fit must contain at least three selectable variables.")
#     }

#     ui <- fluidPage(
#         titlePanel("Interactive Stan posterior explorer"),

#         sidebarLayout(
#             sidebarPanel(
#                 parameterInput(
#                     "x",
#                     "X parameter"
#                 ),
#                 parameterInput(
#                     "y",
#                     "Y parameter"
#                 ),
#                 parameterInput(
#                     "z",
#                     "Z parameter"
#                 ),

#                 radioButtons(
#                     "colour",
#                     "Colour points by",
#                     choices = c("Divergence", "Chain"),
#                     selected = "Divergence"
#                 ),

#                 sliderInput(
#                     "opacity",
#                     "Point opacity",
#                     min = 0.05,
#                     max = 1,
#                     value = 0.6
#                 ),

#                 sliderInput(
#                     "size",
#                     "Point size",
#                     min = 1,
#                     max = 8,
#                     value = 2
#                 )
#             ),

#             mainPanel(
#                 plotlyOutput("posterior_plot", height = "750px")
#             )
#         )
#     )

#     server <- function(input, output, session) {
#         session$onFlushed(
#             function() {
#                 updateSelectizeInput(
#                     session,
#                     "x",
#                     choices = parameters,
#                     selected = parameters[1],
#                     server = TRUE
#                 )
#                 updateSelectizeInput(
#                     session,
#                     "y",
#                     choices = parameters,
#                     selected = parameters[2],
#                     server = TRUE
#                 )
#                 updateSelectizeInput(
#                     session,
#                     "z",
#                     choices = parameters,
#                     selected = parameters[3],
#                     server = TRUE
#                 )
#             },
#             once = TRUE
#         )

#         output$posterior_plot <- renderPlotly({
#             req(input$x, input$y, input$z)

#             colour_variable <- if (input$colour == "Divergence") {
#                 plot_data$divergence
#             } else {
#                 factor(plot_data$.chain)
#             }

#             plot_ly(
#                 x = plot_data[[input$x]],
#                 y = plot_data[[input$y]],
#                 z = plot_data[[input$z]],
#                 color = colour_variable,
#                 type = "scatter3d",
#                 mode = "markers",
#                 marker = list(
#                     size = input$size,
#                     opacity = input$opacity
#                 )
#             ) |>
#                 layout(
#                     scene = list(
#                         xaxis = list(title = input$x),
#                         yaxis = list(title = input$y),
#                         zaxis = list(title = input$z)
#                     ),
#                     legend = list(title = list(text = input$colour))
#                 )
#         })
#     }

#     shinyApp(ui, server)
# }

# launch_fit_comparison <- function(fit1, fit2) {
#     library(shiny)
#     library(plotly)

#     draws1 <- as.data.frame(fit1$draws(format = "df"))
#     draws2 <- as.data.frame(fit2$draws(format = "df"))

#     excluded <- c(".chain", ".iteration", ".draw", "lp__")

#     variables1 <- names(draws1)[
#         !names(draws1) %in% excluded &
#             !grepl("__$", names(draws1))
#     ]

#     variables2 <- names(draws2)[
#         !names(draws2) %in% excluded &
#             !grepl("__$", names(draws2))
#     ]

#     if (length(variables1) < 3 || length(variables2) < 3) {
#         stop("Each fit must contain at least three selectable variables.")
#     }

#     ui <- fluidPage(
#         titlePanel("Compare two Stan posterior distributions"),

#         sidebarLayout(
#             sidebarPanel(
#                 h4("X axis"),
#                 fluidRow(
#                     column(
#                         6,
#                         parameterInput(
#                             "x1",
#                             "Fit 1"
#                         )
#                     ),
#                     column(
#                         6,
#                         parameterInput(
#                             "x2",
#                             "Fit 2"
#                         )
#                     )
#                 ),

#                 h4("Y axis"),
#                 fluidRow(
#                     column(
#                         6,
#                         parameterInput(
#                             "y1",
#                             "Fit 1"
#                         )
#                     ),
#                     column(
#                         6,
#                         parameterInput(
#                             "y2",
#                             "Fit 2"
#                         )
#                     )
#                 ),

#                 h4("Z axis"),
#                 fluidRow(
#                     column(
#                         6,
#                         parameterInput(
#                             "z1",
#                             "Fit 1"
#                         )
#                     ),
#                     column(
#                         6,
#                         parameterInput(
#                             "z2",
#                             "Fit 2"
#                         )
#                     )
#                 ),

#                 sliderInput(
#                     "point_size",
#                     "Point size",
#                     min = 1,
#                     max = 8,
#                     value = 2
#                 ),

#                 sliderInput(
#                     "opacity",
#                     "Opacity",
#                     min = 0.05,
#                     max = 1,
#                     value = 0.45
#                 ),

#                 checkboxInput("show_fit1", "Show Fit 1", value = TRUE),
#                 checkboxInput("show_fit2", "Show Fit 2", value = TRUE)
#             ),

#             mainPanel(
#                 plotlyOutput("posterior_plot", height = "800px")
#             )
#         )
#     )

#     server <- function(input, output, session) {
#         session$onFlushed(
#             function() {
#                 initial_values1 <- variables1[1:3]
#                 initial_values2 <- variables2[1:3]

#                 for (axis_index in seq_along(c("x", "y", "z"))) {
#                     axis <- c("x", "y", "z")[axis_index]

#                     updateSelectizeInput(
#                         session,
#                         paste0(axis, "1"),
#                         choices = variables1,
#                         selected = initial_values1[axis_index],
#                         server = TRUE
#                     )
#                     updateSelectizeInput(
#                         session,
#                         paste0(axis, "2"),
#                         choices = variables2,
#                         selected = initial_values2[axis_index],
#                         server = TRUE
#                     )
#                 }
#             },
#             once = TRUE
#         )

#         observeEvent(input$x1, {
#             if (!is.null(input$x1) && input$x1 %in% variables2) {
#                 updateSelectizeInput(
#                     session,
#                     inputId = "x2",
#                     choices = variables2,
#                     selected = input$x1,
#                     server = TRUE
#                 )
#             }
#         })

#         observeEvent(input$y1, {
#             if (!is.null(input$y1) && input$y1 %in% variables2) {
#                 updateSelectizeInput(
#                     session,
#                     inputId = "y2",
#                     choices = variables2,
#                     selected = input$y1,
#                     server = TRUE
#                 )
#             }
#         })

#         observeEvent(input$z1, {
#             if (!is.null(input$z1) && input$z1 %in% variables2) {
#                 updateSelectizeInput(
#                     session,
#                     inputId = "z2",
#                     choices = variables2,
#                     selected = input$z1,
#                     server = TRUE
#                 )
#             }
#         })

#         output$posterior_plot <- renderPlotly({
#             req(input$x1, input$x2, input$y1, input$y2, input$z1, input$z2)

#             plot <- plot_ly()

#             if (input$show_fit1) {
#                 plot <- plot |>
#                     add_trace(
#                         x = draws1[[input$x1]],
#                         y = draws1[[input$y1]],
#                         z = draws1[[input$z1]],
#                         type = "scatter3d",
#                         mode = "markers",
#                         name = "Fit 1",
#                         marker = list(
#                             size = input$point_size,
#                             opacity = input$opacity,
#                             color = "#2563EB"
#                         ),
#                         text = paste0(
#                             input$x1,
#                             ": ",
#                             signif(draws1[[input$x1]], 4),
#                             "<br>",
#                             input$y1,
#                             ": ",
#                             signif(draws1[[input$y1]], 4),
#                             "<br>",
#                             input$z1,
#                             ": ",
#                             signif(draws1[[input$z1]], 4),
#                             "<br>Chain: ",
#                             draws1$.chain
#                         ),
#                         hoverinfo = "text"
#                     )
#             }

#             if (input$show_fit2) {
#                 plot <- plot |>
#                     add_trace(
#                         x = draws2[[input$x2]],
#                         y = draws2[[input$y2]],
#                         z = draws2[[input$z2]],
#                         type = "scatter3d",
#                         mode = "markers",
#                         name = "Fit 2",
#                         marker = list(
#                             size = input$point_size,
#                             opacity = input$opacity,
#                             color = "#DC2626"
#                         ),
#                         text = paste0(
#                             input$x2,
#                             ": ",
#                             signif(draws2[[input$x2]], 4),
#                             "<br>",
#                             input$y2,
#                             ": ",
#                             signif(draws2[[input$y2]], 4),
#                             "<br>",
#                             input$z2,
#                             ": ",
#                             signif(draws2[[input$z2]], 4),
#                             "<br>Chain: ",
#                             draws2$.chain
#                         ),
#                         hoverinfo = "text"
#                     )
#             }

#             plot |>
#                 layout(
#                     scene = list(
#                         xaxis = list(title = paste(input$x1, "↔", input$x2)),
#                         yaxis = list(title = paste(input$y1, "↔", input$y2)),
#                         zaxis = list(title = paste(input$z1, "↔", input$z2))
#                     ),
#                     legend = list(title = list(text = "Posterior"))
#                 )
#         })
#     }

#     shinyApp(ui, server)
# }

# # parameterInput <- function(inputId, label, choices, selected = NULL) {
# #     shiny::selectizeInput(
# #         inputId = inputId,
# #         label = label,
# #         choices = choices,
# #         selected = selected,
# #         multiple = FALSE,
# #         options = list(
# #             placeholder = "Type to search parameters...",
# #             closeAfterSelect = TRUE,
# #             openOnFocus = TRUE
# #         )
# #     )
# # }

# # launch_stan_3d <- function(fit) {
# #     library(shiny)
# #     library(plotly)
# #     library(dplyr)

# #     draws <- fit$draws(format = "df")

# #     diagnostics <- fit$sampler_diagnostics(format = "df") |>
# #         select(.chain, .iteration, divergent__)

# #     plot_data <- draws |>
# #         left_join(
# #             diagnostics,
# #             by = c(".chain", ".iteration")
# #         ) |>
# #         mutate(
# #             divergence = factor(
# #                 divergent__,
# #                 levels = c(0, 1),
# #                 labels = c("Regular", "Divergent")
# #             )
# #         )

# #     excluded <- c(
# #         ".chain",
# #         ".iteration",
# #         ".draw",
# #         "lp__",
# #         "divergent__",
# #         "divergence"
# #     )

# #     parameters <- names(plot_data)[
# #         !names(plot_data) %in% excluded &
# #             !grepl("__$", names(plot_data))
# #     ]

# #     if (length(parameters) < 3) {
# #         stop("The fit must contain at least three selectable variables.")
# #     }

# #     ui <- fluidPage(
# #         titlePanel("Interactive Stan posterior explorer"),

# #         sidebarLayout(
# #             sidebarPanel(
# #                 parameterInput(
# #                     "x",
# #                     "X parameter",
# #                     parameters,
# #                     selected = parameters[1]
# #                 ),
# #                 parameterInput(
# #                     "y",
# #                     "Y parameter",
# #                     parameters,
# #                     selected = parameters[2]
# #                 ),
# #                 parameterInput(
# #                     "z",
# #                     "Z parameter",
# #                     parameters,
# #                     selected = parameters[3]
# #                 ),

# #                 radioButtons(
# #                     "colour",
# #                     "Colour points by",
# #                     choices = c("Divergence", "Chain"),
# #                     selected = "Divergence"
# #                 ),

# #                 sliderInput(
# #                     "opacity",
# #                     "Point opacity",
# #                     min = 0.05,
# #                     max = 1,
# #                     value = 0.6
# #                 ),

# #                 sliderInput(
# #                     "size",
# #                     "Point size",
# #                     min = 1,
# #                     max = 8,
# #                     value = 2
# #                 )
# #             ),

# #             mainPanel(
# #                 plotlyOutput("posterior_plot", height = "750px")
# #             )
# #         )
# #     )

# #     server <- function(input, output, session) {
# #         output$posterior_plot <- renderPlotly({
# #             req(input$x, input$y, input$z)

# #             colour_variable <- if (input$colour == "Divergence") {
# #                 plot_data$divergence
# #             } else {
# #                 factor(plot_data$.chain)
# #             }

# #             plot_ly(
# #                 x = plot_data[[input$x]],
# #                 y = plot_data[[input$y]],
# #                 z = plot_data[[input$z]],
# #                 color = colour_variable,
# #                 type = "scatter3d",
# #                 mode = "markers",
# #                 marker = list(
# #                     size = input$size,
# #                     opacity = input$opacity
# #                 )
# #             ) |>
# #                 layout(
# #                     scene = list(
# #                         xaxis = list(title = input$x),
# #                         yaxis = list(title = input$y),
# #                         zaxis = list(title = input$z)
# #                     ),
# #                     legend = list(title = list(text = input$colour))
# #                 )
# #         })
# #     }

# #     shinyApp(ui, server)
# # }

# # launch_fit_comparison <- function(fit1, fit2) {
# #     library(shiny)
# #     library(plotly)

# #     draws1 <- fit1$draws(format = "df")
# #     draws2 <- fit2$draws(format = "df")

# #     excluded <- c(".chain", ".iteration", ".draw", "lp__")

# #     variables1 <- names(draws1)[
# #         !names(draws1) %in% excluded &
# #             !grepl("__$", names(draws1))
# #     ]

# #     variables2 <- names(draws2)[
# #         !names(draws2) %in% excluded &
# #             !grepl("__$", names(draws2))
# #     ]

# #     if (length(variables1) < 3 || length(variables2) < 3) {
# #         stop("Each fit must contain at least three selectable variables.")
# #     }

# #     ui <- fluidPage(
# #         titlePanel("Compare two Stan posterior distributions"),

# #         sidebarLayout(
# #             sidebarPanel(
# #                 h4("X axis"),
# #                 fluidRow(
# #                     column(
# #                         6,
# #                         parameterInput(
# #                             "x1",
# #                             "Fit 1",
# #                             variables1,
# #                             selected = variables1[1]
# #                         )
# #                     ),
# #                     column(
# #                         6,
# #                         parameterInput(
# #                             "x2",
# #                             "Fit 2",
# #                             variables2,
# #                             selected = variables2[1]
# #                         )
# #                     )
# #                 ),

# #                 h4("Y axis"),
# #                 fluidRow(
# #                     column(
# #                         6,
# #                         parameterInput(
# #                             "y1",
# #                             "Fit 1",
# #                             variables1,
# #                             selected = variables1[2]
# #                         )
# #                     ),
# #                     column(
# #                         6,
# #                         parameterInput(
# #                             "y2",
# #                             "Fit 2",
# #                             variables2,
# #                             selected = variables2[2]
# #                         )
# #                     )
# #                 ),

# #                 h4("Z axis"),
# #                 fluidRow(
# #                     column(
# #                         6,
# #                         parameterInput(
# #                             "z1",
# #                             "Fit 1",
# #                             variables1,
# #                             selected = variables1[3]
# #                         )
# #                     ),
# #                     column(
# #                         6,
# #                         parameterInput(
# #                             "z2",
# #                             "Fit 2",
# #                             variables2,
# #                             selected = variables2[3]
# #                         )
# #                     )
# #                 ),

# #                 sliderInput(
# #                     "point_size",
# #                     "Point size",
# #                     min = 1,
# #                     max = 8,
# #                     value = 2
# #                 ),

# #                 sliderInput(
# #                     "opacity",
# #                     "Opacity",
# #                     min = 0.05,
# #                     max = 1,
# #                     value = 0.45
# #                 ),

# #                 checkboxInput("show_fit1", "Show Fit 1", value = TRUE),
# #                 checkboxInput("show_fit2", "Show Fit 2", value = TRUE)
# #             ),

# #             mainPanel(
# #                 plotlyOutput("posterior_plot", height = "800px")
# #             )
# #         )
# #     )

# #     server <- function(input, output, session) {
# #         observeEvent(input$x1, {
# #             if (!is.null(input$x1) && input$x1 %in% variables2) {
# #                 updateSelectizeInput(
# #                     session,
# #                     inputId = "x2",
# #                     selected = input$x1
# #                 )
# #             }
# #         })

# #         observeEvent(input$y1, {
# #             if (!is.null(input$y1) && input$y1 %in% variables2) {
# #                 updateSelectizeInput(
# #                     session,
# #                     inputId = "y2",
# #                     selected = input$y1
# #                 )
# #             }
# #         })

# #         observeEvent(input$z1, {
# #             if (!is.null(input$z1) && input$z1 %in% variables2) {
# #                 updateSelectizeInput(
# #                     session,
# #                     inputId = "z2",
# #                     selected = input$z1
# #                 )
# #             }
# #         })

# #         output$posterior_plot <- renderPlotly({
# #             req(input$x1, input$x2, input$y1, input$y2, input$z1, input$z2)

# #             plot <- plot_ly()

# #             if (input$show_fit1) {
# #                 plot <- plot |>
# #                     add_trace(
# #                         x = draws1[[input$x1]],
# #                         y = draws1[[input$y1]],
# #                         z = draws1[[input$z1]],
# #                         type = "scatter3d",
# #                         mode = "markers",
# #                         name = "Fit 1",
# #                         marker = list(
# #                             size = input$point_size,
# #                             opacity = input$opacity,
# #                             color = "#2563EB"
# #                         ),
# #                         text = paste0(
# #                             input$x1,
# #                             ": ",
# #                             signif(draws1[[input$x1]], 4),
# #                             "<br>",
# #                             input$y1,
# #                             ": ",
# #                             signif(draws1[[input$y1]], 4),
# #                             "<br>",
# #                             input$z1,
# #                             ": ",
# #                             signif(draws1[[input$z1]], 4),
# #                             "<br>Chain: ",
# #                             draws1$.chain
# #                         ),
# #                         hoverinfo = "text"
# #                     )
# #             }

# #             if (input$show_fit2) {
# #                 plot <- plot |>
# #                     add_trace(
# #                         x = draws2[[input$x2]],
# #                         y = draws2[[input$y2]],
# #                         z = draws2[[input$z2]],
# #                         type = "scatter3d",
# #                         mode = "markers",
# #                         name = "Fit 2",
# #                         marker = list(
# #                             size = input$point_size,
# #                             opacity = input$opacity,
# #                             color = "#DC2626"
# #                         ),
# #                         text = paste0(
# #                             input$x2,
# #                             ": ",
# #                             signif(draws2[[input$x2]], 4),
# #                             "<br>",
# #                             input$y2,
# #                             ": ",
# #                             signif(draws2[[input$y2]], 4),
# #                             "<br>",
# #                             input$z2,
# #                             ": ",
# #                             signif(draws2[[input$z2]], 4),
# #                             "<br>Chain: ",
# #                             draws2$.chain
# #                         ),
# #                         hoverinfo = "text"
# #                     )
# #             }

# #             plot |>
# #                 layout(
# #                     scene = list(
# #                         xaxis = list(title = paste(input$x1, "↔", input$x2)),
# #                         yaxis = list(title = paste(input$y1, "↔", input$y2)),
# #                         zaxis = list(title = paste(input$z1, "↔", input$z2))
# #                     ),
# #                     legend = list(title = list(text = "Posterior"))
# #                 )
# #         })
# #     }

# #     shinyApp(ui, server)
# # }
