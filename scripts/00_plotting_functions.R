# ---- tabular-io ----
read_tsv_table <- function(file) {
  utils::read.delim(
    file,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_tsv_table <- function(x, file) {
  utils::write.table(
    x,
    file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# ---- device-helpers ----
with_pdf_device <- function(file, width, height, expr) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = width, height = height)
  old_par <- graphics::par(no.readonly = TRUE)

  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  force(expr)
  invisible(file)
}

# ---- layout-helpers ----
model_panel_layout <- function(n, max_cols = 4) {
  if (n < 1) {
    stop("n must be at least 1.")
  }

  ncol <- min(max_cols, n)
  nrow <- ceiling(n / ncol)
  c(nrow = nrow, ncol = ncol)
}

model_panel_pdf_height <- function(n, row_height = 2.4, min_height = 7, max_cols = 4) {
  layout_dim <- model_panel_layout(n, max_cols = max_cols)
  max(min_height, row_height * layout_dim["nrow"])
}

facet_panel_layout <- function(n, max_cols = 4) {
  if (n <= 1) {
    return(c(nrow = 1, ncol = 1))
  }

  if (n == 2) {
    return(c(nrow = 1, ncol = 2))
  }

  if (n <= 4) {
    return(c(nrow = 2, ncol = 2))
  }

  ncol <- min(max_cols, ceiling(sqrt(n)))
  nrow <- ceiling(n / ncol)
  c(nrow = nrow, ncol = ncol)
}

facet_panel_pdf_size <- function(
  n,
  width = 12,
  one_or_two_height = 4.5,
  four_panel_height = 8,
  row_height = 2.4,
  min_height = 8
) {
  layout_dim <- facet_panel_layout(n)

  if (n <= 2) {
    return(c(width = width, height = one_or_two_height))
  }

  if (n <= 4) {
    return(c(width = width, height = four_panel_height))
  }

  c(width = width, height = max(min_height, row_height * layout_dim["nrow"]))
}

plot_range <- function(x, pad_fraction = 0.04) {
  out <- range(x, na.rm = TRUE)

  if (!all(is.finite(out))) {
    return(c(0, 1))
  }

  if (diff(out) == 0) {
    pad <- if (out[1] == 0) 0.5 else abs(out[1]) * 0.05
    return(out + c(-pad, pad))
  }

  pad <- diff(out) * pad_fraction
  out + c(-pad, pad)
}

file_stub <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", x)
  out <- gsub("^_+|_+$", "", out)
  tolower(out)
}

empty_panel <- function(main = "No data") {
  graphics::plot.new()
  graphics::title(main = main)
}

# ---- style-helpers ----
response_label <- function(response) {
  toupper(sub("[.].*$", "", response))
}

make_response_plot_styles <- function(responses, labels = NULL) {
  if (is.null(labels)) {
    labels <- response_label(responses)
  }

  pch_values <- c(16, 2, 17, 15, 1, 0)
  line_cols <- c("black", "#D95F5F", "#1B9E77", "#7570B3", "#E6AB02", "#666666")
  line_ltys <- c(1, 2, 3, 4, 5, 6)

  data.frame(
    response = responses,
    label = labels,
    pch = pch_values[seq_along(responses)],
    line_col = line_cols[seq_along(responses)],
    line_lty = line_ltys[seq_along(responses)],
    stringsAsFactors = FALSE
  )
}

response_style <- function(styles, response) {
  styles[match(response, styles$response), , drop = FALSE]
}

draw_response_legend <- function(styles, position = "topright") {
  graphics::legend(
    position,
    legend = styles$label,
    col = styles$line_col,
    pch = styles$pch,
    lty = styles$line_lty,
    lwd = 2.4,
    bty = "n",
    cex = 0.85
  )
}

discrete_palette <- function(levels, palette = "Dark 3") {
  levels <- sort(unique(as.character(levels)))

  if (length(levels) == 0) {
    return(setNames(character(0), character(0)))
  }

  setNames(grDevices::hcl.colors(length(levels), palette = palette), levels)
}

max_panel_count <- function(data, models, responses, panel_col = "h1") {
  counts <- vapply(models, function(model) {
    d <- data[data$model == model & data$response %in% responses, , drop = FALSE]

    if (nrow(d) == 0) {
      return(1L)
    }

    by_response <- split(d, d$response, drop = TRUE)
    max(vapply(by_response, function(x) length(unique(x[[panel_col]])), integer(1)))
  }, integer(1))

  max(counts, na.rm = TRUE)
}

# ---- model-diagnostic-plots ----
plot_observed_vs_predicted_pdf <- function(
  data,
  file,
  responses = unique(data$response),
  models = unique(data$model),
  width = 12
) {
  layout_dim <- model_panel_layout(length(models))

  with_pdf_device(file, width = width, height = model_panel_pdf_height(length(models)), {
    for (resp in responses) {
      graphics::par(
        mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
        mar = c(4, 4, 3, 1),
        oma = c(0, 0, 2, 0)
      )

      for (model in models) {
        d <- data[data$response == resp & data$model == model, , drop = FALSE]

        if (nrow(d) == 0) {
          empty_panel(paste("No data:", resp, model))
          next
        }

        lim <- plot_range(c(d$observed, d$predicted))

        graphics::plot(
          d$observed,
          d$predicted,
          xlim = lim,
          ylim = lim,
          pch = 16,
          col = grDevices::rgb(0, 0, 0, 0.35),
          xlab = "Observed",
          ylab = "Predicted",
          main = paste(resp, model, sep = " | ")
        )

        graphics::abline(0, 1, col = "red", lwd = 2)
      }

      graphics::mtext(paste("Observed vs predicted:", resp), outer = TRUE, cex = 1.2)
    }
  })
}

plot_response_vs_x_pdf <- function(
  observed,
  predicted,
  file,
  responses = unique(observed$response),
  models = unique(observed$model),
  x_col = "x",
  x_lab = "Relative stand density (rsd)",
  width = 12
) {
  layout_dim <- model_panel_layout(length(models))

  with_pdf_device(file, width = width, height = model_panel_pdf_height(length(models)), {
    for (resp in responses) {
      graphics::par(
        mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
        mar = c(4, 4, 3, 1),
        oma = c(0, 0, 2, 0)
      )

      for (model in models) {
        obs_d <- observed[observed$response == resp & observed$model == model, , drop = FALSE]
        cur_d <- predicted[predicted$response == resp & predicted$model == model, , drop = FALSE]

        if (nrow(obs_d) == 0 || nrow(cur_d) == 0) {
          empty_panel(paste("No data:", resp, model))
          next
        }

        ord <- order(cur_d[[x_col]])
        ylim <- plot_range(c(obs_d$observed, cur_d$pred_q05, cur_d$pred_q95))

        graphics::plot(
          obs_d[[x_col]],
          obs_d$observed,
          pch = 16,
          col = grDevices::rgb(0, 0, 0, 0.35),
          xlab = x_lab,
          ylab = resp,
          main = paste(resp, model, sep = " | "),
          ylim = ylim
        )

        graphics::lines(cur_d[[x_col]][ord], cur_d$predicted[ord], col = "red", lwd = 2)
        graphics::lines(cur_d[[x_col]][ord], cur_d$pred_q05[ord], col = "red", lty = 2)
        graphics::lines(cur_d[[x_col]][ord], cur_d$pred_q95[ord], col = "red", lty = 2)
      }

      graphics::mtext(paste("Response vs rsd:", resp), outer = TRUE, cex = 1.2)
    }
  })
}

# ---- grouped-response-plots ----
draw_grouped_response_page <- function(
  model,
  observed,
  predicted,
  responses,
  styles,
  panel_col = "h1",
  color_col = "h2",
  x_col = "x",
  x_lab = "Relative stand density",
  y_lab = "Biomass expansion factors"
) {
  obs_d <- observed[observed$model == model & observed$response %in% responses, , drop = FALSE]
  cur_d <- predicted[predicted$model == model & predicted$response %in% responses, , drop = FALSE]

  if (nrow(obs_d) == 0 || nrow(cur_d) == 0) {
    empty_panel(paste("No data:", model))
    return(invisible(NULL))
  }

  panel_levels <- sort(unique(obs_d[[panel_col]]))
  color_levels <- sort(unique(c(obs_d[[color_col]], cur_d[[color_col]])))
  color_map <- discrete_palette(color_levels)
  layout_dim <- facet_panel_layout(length(panel_levels))
  layout_slots <- prod(layout_dim)
  xlim <- plot_range(c(obs_d[[x_col]], cur_d[[x_col]]))
  ylim <- plot_range(c(obs_d$observed, cur_d$predicted))

  graphics::par(
    mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (panel_i in seq_along(panel_levels)) {
    panel_level <- panel_levels[[panel_i]]

    graphics::plot(
      NA,
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = x_lab,
      ylab = y_lab,
      main = panel_level
    )

    for (resp in responses) {
      style <- response_style(styles, resp)

      for (color_level in color_levels) {
        obs_i <- obs_d[
          obs_d$response == resp & obs_d[[panel_col]] == panel_level & obs_d[[color_col]] == color_level,
          ,
          drop = FALSE
        ]
        cur_i <- cur_d[
          cur_d$response == resp & cur_d[[panel_col]] == panel_level & cur_d[[color_col]] == color_level,
          ,
          drop = FALSE
        ]
        col_i <- color_map[[color_level]]

        if (nrow(obs_i) > 0) {
          graphics::points(
            obs_i[[x_col]],
            obs_i$observed,
            pch = style$pch,
            col = grDevices::adjustcolor(col_i, alpha.f = 0.55)
          )
        }

        if (nrow(cur_i) > 0) {
          ord <- order(cur_i[[x_col]])
          graphics::lines(
            cur_i[[x_col]][ord],
            cur_i$predicted[ord],
            col = grDevices::adjustcolor(style$line_col, alpha.f = 0.80),
            lwd = 2.4,
            lty = style$line_lty
          )
        }
      }
    }

    if (panel_i == 1) {
      draw_response_legend(styles)
    }
  }

  if (layout_slots > length(panel_levels)) {
    for (i in seq_len(layout_slots - length(panel_levels))) {
      graphics::plot.new()
    }
  }

  graphics::mtext(paste("Response vs rsd by h1/h2:", model), outer = TRUE, cex = 1.2)
  invisible(NULL)
}

plot_grouped_response_pdf <- function(
  observed,
  predicted,
  file,
  models,
  responses,
  styles = make_response_plot_styles(responses)
) {
  pdf_size <- facet_panel_pdf_size(max_panel_count(observed, models, responses))

  with_pdf_device(file, width = pdf_size["width"], height = pdf_size["height"], {
    for (model in models) {
      draw_grouped_response_page(model, observed, predicted, responses, styles)
    }
  })
}

plot_grouped_response_pdfs_by_model <- function(
  observed,
  predicted,
  out_dir,
  file_prefix,
  models,
  responses,
  styles = make_response_plot_styles(responses)
) {
  vapply(models, function(model) {
    n_panels <- max_panel_count(observed, model, responses)
    pdf_size <- facet_panel_pdf_size(n_panels)
    file <- paste0(file_prefix, file_stub(model), ".pdf")

    with_pdf_device(file.path(out_dir, file), width = pdf_size["width"], height = pdf_size["height"], {
      draw_grouped_response_page(model, observed, predicted, responses, styles)
    })

    file
  }, character(1))
}

# ---- grouped-residual-plots ----
draw_grouped_residual_page <- function(
  model,
  data,
  responses,
  styles,
  x_var,
  x_lab,
  panel_col = "h1",
  color_col = "h2"
) {
  d <- data[data$model == model & data$response %in% responses, , drop = FALSE]

  if (nrow(d) == 0) {
    empty_panel(paste("No residuals:", model))
    return(invisible(NULL))
  }

  panel_levels <- sort(unique(d[[panel_col]]))
  color_levels <- sort(unique(d[[color_col]]))
  color_map <- discrete_palette(color_levels)
  layout_dim <- facet_panel_layout(length(panel_levels))
  layout_slots <- prod(layout_dim)
  xlim <- plot_range(d[[x_var]])
  ylim <- plot_range(d$residual)

  graphics::par(
    mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (panel_i in seq_along(panel_levels)) {
    panel_level <- panel_levels[[panel_i]]

    graphics::plot(
      NA,
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = x_lab,
      ylab = "Residual",
      main = panel_level
    )
    graphics::abline(h = 0, col = "gray45", lty = 2)

    for (resp in responses) {
      style <- response_style(styles, resp)

      for (color_level in color_levels) {
        d_i <- d[
          d$response == resp & d[[panel_col]] == panel_level & d[[color_col]] == color_level,
          ,
          drop = FALSE
        ]
        col_i <- color_map[[color_level]]

        if (nrow(d_i) > 0) {
          graphics::points(
            d_i[[x_var]],
            d_i$residual,
            pch = style$pch,
            col = grDevices::adjustcolor(col_i, alpha.f = 0.55)
          )
        }
      }
    }

    if (panel_i == 1) {
      draw_response_legend(styles)
    }
  }

  if (layout_slots > length(panel_levels)) {
    for (i in seq_len(layout_slots - length(panel_levels))) {
      graphics::plot.new()
    }
  }

  graphics::mtext(paste("Residual vs", x_lab, "by h1/h2:", model), outer = TRUE, cex = 1.2)
  invisible(NULL)
}

plot_grouped_residual_pdf <- function(
  data,
  file,
  models,
  responses,
  styles,
  x_var,
  x_lab
) {
  pdf_size <- facet_panel_pdf_size(max_panel_count(data, models, responses))

  with_pdf_device(file, width = pdf_size["width"], height = pdf_size["height"], {
    for (model in models) {
      draw_grouped_residual_page(model, data, responses, styles, x_var = x_var, x_lab = x_lab)
    }
  })
}

plot_grouped_residual_pdfs_by_model <- function(
  data,
  out_dir,
  file_prefix,
  models,
  responses,
  styles
) {
  vapply(models, function(model) {
    n_panels <- max_panel_count(data, model, responses)
    pdf_size <- facet_panel_pdf_size(n_panels)
    file <- paste0(file_prefix, file_stub(model), ".pdf")

    with_pdf_device(file.path(out_dir, file), width = pdf_size["width"], height = pdf_size["height"], {
      draw_grouped_residual_page(model, data, responses, styles, x_var = "predicted", x_lab = "Predicted")
      draw_grouped_residual_page(model, data, responses, styles, x_var = "x", x_lab = "Relative stand density")
    })

    file
  }, character(1))
}

# ---- grouped-mean-residual-plots ----
draw_grouped_mean_residual_page <- function(
  model,
  data,
  responses,
  styles,
  mean_group,
  panel_col = "h1",
  color_col = "h2"
) {
  d <- data[
    data$model == model &
      data$response %in% responses &
      data$mean_group == mean_group,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0) {
    empty_panel(paste("No mean residuals:", mean_group, "|", model))
    return(invisible(NULL))
  }

  panel_levels <- sort(unique(d[[panel_col]]))
  color_levels <- sort(unique(d[[color_col]]))
  color_map <- discrete_palette(color_levels)
  layout_dim <- facet_panel_layout(length(panel_levels))
  layout_slots <- prod(layout_dim)
  xlim <- plot_range(d$x)
  ylim <- plot_range(d$residual_from_mean)

  graphics::par(
    mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (panel_i in seq_along(panel_levels)) {
    panel_level <- panel_levels[[panel_i]]

    graphics::plot(
      NA,
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = "Relative stand density",
      ylab = paste("Residual from", mean_group, "mean"),
      main = panel_level
    )
    graphics::abline(h = 0, col = "gray45", lty = 2)

    for (resp in responses) {
      style <- response_style(styles, resp)

      for (color_level in color_levels) {
        d_i <- d[
          d$response == resp & d[[panel_col]] == panel_level & d[[color_col]] == color_level,
          ,
          drop = FALSE
        ]
        col_i <- color_map[[color_level]]

        if (nrow(d_i) > 0) {
          graphics::points(
            d_i$x,
            d_i$residual_from_mean,
            pch = style$pch,
            col = grDevices::adjustcolor(col_i, alpha.f = 0.55)
          )
        }
      }
    }

    if (panel_i == 1) {
      draw_response_legend(styles)
    }
  }

  if (layout_slots > length(panel_levels)) {
    for (i in seq_len(layout_slots - length(panel_levels))) {
      graphics::plot.new()
    }
  }

  graphics::mtext(
    paste("Residual from", mean_group, "mean vs rsd by h1/h2:", model),
    outer = TRUE,
    cex = 1.2
  )
  invisible(NULL)
}

plot_grouped_mean_residual_pdf <- function(
  data,
  file,
  models,
  responses,
  styles,
  mean_group
) {
  pdf_size <- facet_panel_pdf_size(max_panel_count(data, models, responses))

  with_pdf_device(file, width = pdf_size["width"], height = pdf_size["height"], {
    for (model in models) {
      draw_grouped_mean_residual_page(model, data, responses, styles, mean_group = mean_group)
    }
  })
}

plot_grouped_mean_residual_pdfs_by_model <- function(
  data,
  out_dir,
  file_prefix,
  models,
  responses,
  styles,
  mean_groups = unique(data$mean_group)
) {
  vapply(models, function(model) {
    n_panels <- max_panel_count(data, model, responses)
    pdf_size <- facet_panel_pdf_size(n_panels)
    file <- paste0(file_prefix, file_stub(model), ".pdf")

    with_pdf_device(file.path(out_dir, file), width = pdf_size["width"], height = pdf_size["height"], {
      for (mean_group in mean_groups) {
        draw_grouped_mean_residual_page(model, data, responses, styles, mean_group = mean_group)
      }
    })

    file
  }, character(1))
}
