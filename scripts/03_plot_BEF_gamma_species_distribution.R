# Species-level observed BEF distributions for the mono_N / mono_B comparison.
# Points show observations; the solid curves are method-of-moments Gamma fits.

current_script_dir <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files)]

  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[[length(frame_files)]], mustWork = TRUE)))
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }

  getwd()
}

script_dir <- current_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
setwd(project_root)
source(file.path(script_dir, "00_plotting_functions.R"))

input_file <- file.path(project_root, "processed_data", "plot_biomass.txt")
figure_dir <- file.path(project_root, "figures", "Bayes_gamma")
output_dir <- file.path(project_root, "processed_data", "bef_bayes_gamma_species")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

bp <- read_tsv_table(input_file)
required_cols <- c(
  "sp_code", "species", "Family", "Genus", "ft1.forest_type", "befa.st", "befr.st"
)
missing_cols <- setdiff(required_cols, names(bp))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

forest_types <- c("mono_N", "mono_B")
responses <- c("befa.st", "befr.st")
response_labels <- c("befa.st" = "BEFA", "befr.st" = "BEFR")

bp <- bp[bp$ft1.forest_type %in% forest_types, , drop = FALSE]
bp$sp_code <- as.character(bp$sp_code)
bp$species <- as.character(bp$species)

species_metadata <- unique(bp[, c(
  "sp_code", "species", "ft1.forest_type", "Family", "Genus"
)])
species_order <- split(species_metadata, species_metadata$ft1.forest_type)
species_order <- lapply(species_order, function(x) {
  x <- x[order(x$Family, x$Genus, x$species, x$sp_code), , drop = FALSE]
  x$sp_code
})

species_metadata$species_key <- paste(species_metadata$ft1.forest_type, species_metadata$sp_code, sep = "|")
species_colors <- lapply(forest_types, function(forest_type) {
  codes <- species_order[[forest_type]]
  setNames(grDevices::rainbow(length(codes), start = 0.02, end = 0.94), codes)
})
names(species_colors) <- forest_types

fit_gamma <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2 || stats::var(x) <= 0) {
    return(c(n = length(x), mean = mean(x), variance = stats::var(x), shape = NA, rate = NA))
  }

  x_mean <- mean(x)
  x_var <- stats::var(x)
  c(
    n = length(x),
    mean = x_mean,
    variance = x_var,
    shape = x_mean^2 / x_var,
    rate = x_mean / x_var
  )
}

gamma_rows <- list()
row_i <- 1L
for (forest_type in forest_types) {
  for (sp_code in species_order[[forest_type]]) {
    for (response in responses) {
      x <- bp[bp$ft1.forest_type == forest_type & bp$sp_code == sp_code, response]
      pars <- fit_gamma(x)
      species_name <- unique(bp$species[bp$ft1.forest_type == forest_type & bp$sp_code == sp_code])
      gamma_rows[[row_i]] <- data.frame(
        forest_type = forest_type,
        sp_code = sp_code,
        species = if (length(species_name) == 0) NA_character_ else species_name[[1]],
        response = response,
        n = unname(pars[["n"]]),
        mean = unname(pars[["mean"]]),
        variance = unname(pars[["variance"]]),
        shape = unname(pars[["shape"]]),
        rate = unname(pars[["rate"]]),
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
}

gamma_parameters <- do.call(rbind, gamma_rows)
write_tsv_table(
  gamma_parameters,
  file.path(output_dir, "01_gamma_species_parameters.txt")
)

plot_distribution_panel <- function(
  data,
  forest_type,
  response,
  y_limits,
  panel_label,
  show_x_labels = TRUE
) {
  species <- species_order[[forest_type]]
  panel <- data[data$forest_type == forest_type & data$response == response, , drop = FALSE]
  panel_bp <- bp[
    bp$ft1.forest_type == forest_type &
      is.finite(bp[[response]]) &
      bp$sp_code %in% species,
    , drop = FALSE
  ]
  panel_bp$sp_code <- factor(panel_bp$sp_code, levels = species)
  panel_bp$x_position <- match(as.character(panel_bp$sp_code), species)
  x_positions <- seq_along(species)
  point_pch <- if (response == "befa.st") {
    if (forest_type == "mono_N") 17 else 16
  } else {
    if (forest_type == "mono_N") 25 else 15
  }

  graphics::plot(
    panel_bp$x_position,
    panel_bp[[response]],
    type = "n",
    xlim = c(0, length(species) + 1),
    ylim = y_limits,
    xaxt = "n",
    xlab = "",
    ylab = "",
    yaxt = "n",
    bty = "n",
    lwd = 0.3,
    col = 1
  )

  point_colors <- unname(species_colors[[forest_type]][as.character(panel_bp$sp_code)])
  graphics::points(
    jitter(panel_bp$x_position, amount = 0.11),
    panel_bp[[response]],
    pch = point_pch,
    col = grDevices::adjustcolor(point_colors, alpha.f = 0.05),
    bg = grDevices::adjustcolor(point_colors, alpha.f = 0.05),
    cex = 0.42,
    lwd = 0.15
  )

  for (i in seq_along(species)) {
    sp <- species[[i]]

    pars <- panel[panel$sp_code == sp, , drop = FALSE]
    if (nrow(pars) == 1 && is.finite(pars$shape) && is.finite(pars$rate)) {
      y_grid <- seq(y_limits[[1]], y_limits[[2]], length.out = 250)
      density_values <- stats::dgamma(y_grid, shape = pars$shape, rate = pars$rate)
      density_values <- density_values / max(density_values, na.rm = TRUE)
      width <- 0.62
      curve_col <- species_colors[[forest_type]][[sp]]
      curve_x_left <- i - width * density_values
      curve_x_right <- i + width * density_values
      graphics::polygon(
        c(curve_x_left, rev(curve_x_right)),
        c(y_grid, rev(y_grid)),
        border = "black",
        col = grDevices::adjustcolor(curve_col, alpha.f = 0.25),
        lwd = 2.5
      )
      graphics::lines(
        curve_x_right,
        y_grid,
        col = "black",
        lwd = 2.8
      )
      graphics::lines(
        curve_x_left,
        y_grid,
        col = "black",
        lwd = 2.8
      )
      graphics::lines(
        curve_x_right,
        y_grid,
        col = curve_col,
        lwd = 1.2
      )
      graphics::lines(
        curve_x_left,
        y_grid,
        col = curve_col,
        lwd = 1.2
      )
    }
  }

  graphics::axis(
    1,
    at = x_positions,
    labels = if (show_x_labels) species else FALSE,
    las = 2,
    cex.axis = 7 / 12,
    lwd = 0.3,
    tck = -0.02,
    mgp = c(0, -0.30, 0),
    font = 1
  )
  graphics::axis(
    2,
    at = pretty(y_limits),
    labels = TRUE,
    las = 1,
    cex.axis = 7 / 12,
    lwd = 0.3,
    tck = -0.02,
    mgp = c(0, 0, 0)
  )
  if (forest_type == "mono_N") {
    graphics::mtext(
      response_labels[[response]],
      side = 2,
      line = 0.50,
      font = 1,
      cex = 7 / 12
    )
  }
  graphics::mtext(
    panel_label,
    side = 3,
    line = -0.80,
    adj = 0.02,
    font = 2,
    cex = 8 / 12
  )
  graphics::mtext(
    forest_type,
    side = 3,
    line = -0.10,
    adj = 0.98,
    font = 1,
    cex = 7 / 12
  )
}

y_limits <- setNames(lapply(responses, function(response) {
  x <- bp[[response]]
  x <- x[is.finite(x)]
  c(0, max(x, na.rm = TRUE) * 1.05)
}), responses)

figure_file <- file.path(
  figure_dir,
  "04_gamma_species_distribution_monoN_monoB.pdf"
)

with_pdf_device(figure_file, width = 3.42, height = 3.60, {
  graphics::par(
    mfrow = c(2, 2),
    mai = c(0.40, 0.40, 0.20, 0.02),
    oma = c(0, 0, 0, 0),
    lwd = 0.3,
    col = 1,
    fg = 1,
    bg = "white"
  )

  plot_distribution_panel(gamma_parameters, "mono_N", "befa.st", y_limits[["befa.st"]], "a")
  plot_distribution_panel(gamma_parameters, "mono_B", "befa.st", y_limits[["befa.st"]], "b")
  plot_distribution_panel(gamma_parameters, "mono_N", "befr.st", y_limits[["befr.st"]], "c")
  plot_distribution_panel(gamma_parameters, "mono_B", "befr.st", y_limits[["befr.st"]], "d")
})

writeLines(
  c(
    "Species-level observed BEF Gamma distribution figure",
    paste("Generated:", as.character(Sys.time())),
    paste("Source:", input_file),
    paste("Figure:", figure_file),
    paste("Parameters:", file.path(output_dir, "01_gamma_species_parameters.txt")),
    "The plotting data are read as bp from processed_data/plot_biomass.txt, matching the Figure_style.R workflow.",
    "Points are observed BEF values with low opacity and species-specific symbols.",
    "Curves are method-of-moments Gamma densities, normalized within species for display.",
    "Species are ordered by Family, Genus, species name, and sp_code.",
    "Panels: columns mono_N and mono_B; rows befa.st and befr.st."
  ),
  file.path(output_dir, "README_FIGURE.txt")
)

figure_file
