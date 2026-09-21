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

        output$posterior_plot <- renderPlotly({
            req(input$x, input$y, input$z)

            colour_variable <- if (input$colour == "Divergence") {
                plot_data$divergence
            } else {
                factor(plot_data$.chain)
            }

            plot_ly(
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
            ) |>
                layout(
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
                checkboxInput("show_fit2", "Show Fit 2", value = TRUE)
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

        output$posterior_plot <- renderPlotly({
            plot <- plot_ly()
            visible_fits <- 0L

            for (i in seq_len(number_of_fits)) {
                if (!isTRUE(input[[paste0("show_", i)]])) {
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

        output$draws_posterior_plot <- plotly::renderPlotly({
            plot_data <- selected_draws()

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

            plotly::plot_ly(
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
            ) |>
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
