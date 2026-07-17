# ---- setup ----
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
project_root <- if (basename(script_dir) == "scripts") {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
} else {
  normalizePath(script_dir, mustWork = TRUE)
}

setwd(project_root)
source(file.path(script_dir, "00_plotting_functions.R"))

# ---- directories ----
gamma_output_dir <- file.path(project_root, "processed_data", "bef_bayes_gamma")
source_data_file <- file.path(project_root, "processed_data", "plot_biomass.txt")
analysis_dir <- file.path(project_root, "processed_data", "bef_bayes_gamma_dominant_species")
figure_dir <- file.path(project_root, "figures", "Bayes_gamma_dominant_species")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

split_var <- Sys.getenv("BEF_DOMINANT_SPLIT", "dominant_genus")
min_split_n <- as.integer(Sys.getenv("BEF_DOMINANT_MIN_N", "1"))

if (is.na(min_split_n) || min_split_n < 1) {
  stop("BEF_DOMINANT_MIN_N must be a positive integer.")
}

read_output <- function(file, dir = gamma_output_dir) {
  path <- file.path(dir, file)

  if (!file.exists(path)) {
    stop("Missing required output: ", path)
  }

  read_tsv_table(path)
}

write_output <- function(x, file) {
  write_tsv_table(x, file.path(analysis_dir, file))
}

# ---- data ----
model_info <- read_output("00_model_info.txt")
step5_observed_fit <- read_output("05_observed_vs_predicted.txt")

if (!file.exists(source_data_file)) {
  stop("Missing source observation data: ", source_data_file)
}

source_data <- read_tsv_table(source_data_file)

required_source_cols <- c(
  "stand_id", "plot", "sp_code", "species", "PFT", "Genus", "Family",
  "species.1", "ft1.dominant_prop", "rsd", "befa.st", "befr.st", "beft.st"
)
missing_source_cols <- setdiff(required_source_cols, names(source_data))

if (length(missing_source_cols) > 0) {
  stop("Missing required source columns: ", paste(missing_source_cols, collapse = ", "))
}

response_source_specs <- list(
  befa.st = list(keep = !is.na(source_data$befa.st) & !is.na(source_data$rsd), observed = source_data$befa.st),
  befr.st = list(keep = !is.na(source_data$befr.st) & !is.na(source_data$rsd), observed = source_data$befr.st),
  beft.st = list(keep = !is.na(source_data$befa.st) & !is.na(source_data$befr.st) & !is.na(source_data$rsd), observed = source_data$beft.st)
)

make_response_metadata <- function(response) {
  spec <- response_source_specs[[response]]
  d <- source_data[spec$keep, , drop = FALSE]
  out <- data.frame(
    response = response,
    row_id = seq_len(nrow(d)),
    source_row_id = which(spec$keep),
    source_observed = spec$observed[spec$keep],
    stand_id = d$stand_id,
    plot = d$plot,
    sp_code = d$sp_code,
    tree_species = d$species,
    tree_genus = d$Genus,
    tree_family = d$Family,
    PFT = d$PFT,
    dominant_species = d$species.1,
    dominant_genus = sub("_.*$", "", d$species.1),
    dominant_prop = d$ft1.dominant_prop,
    rsd_source = d$rsd,
    stringsAsFactors = FALSE
  )
  out$dominant_species[is.na(out$dominant_species) | out$dominant_species == ""] <- "unknown"
  out$dominant_genus[is.na(out$dominant_genus) | out$dominant_genus == ""] <- "unknown"
  out
}

response_metadata <- do.call(rbind, lapply(names(response_source_specs), make_response_metadata))

residual_data <- merge(
  step5_observed_fit,
  response_metadata,
  by = c("response", "row_id"),
  all.x = TRUE,
  sort = FALSE
)

missing_metadata <- residual_data[is.na(residual_data$dominant_genus), c("model", "response", "row_id")]
if (nrow(missing_metadata) > 0) {
  stop("Some residual rows could not be matched back to source observations.")
}

max_observed_difference <- max(abs(residual_data$observed - residual_data$source_observed), na.rm = TRUE)
if (is.finite(max_observed_difference) && max_observed_difference > 1e-8) {
  warning("Largest observed/source observed mismatch: ", signif(max_observed_difference, 4))
}

if (!split_var %in% names(residual_data)) {
  stop(
    "Unknown split variable: ", split_var,
    ". Available split variables include dominant_genus, dominant_species, tree_genus, tree_species, PFT, sp_code."
  )
}

add_split_mean_residuals <- function(d, split_col) {
  split_key <- interaction(d$model, d$response, d[[split_col]], drop = TRUE)
  d$mean_observed_split <- ave(d$observed, split_key, FUN = function(x) mean(x, na.rm = TRUE))
  d$residual_from_split_mean <- d$observed - d$mean_observed_split

  for (mean_group in c("h1", "h2")) {
    group_key <- interaction(
      d$model,
      d$response,
      d[[split_col]],
      d[[mean_group]],
      drop = TRUE
    )
    mean_col <- paste0("mean_observed_", mean_group, "_within_split")
    residual_col <- paste0("residual_from_", mean_group, "_mean_within_split")
    d[[mean_col]] <- ave(d$observed, group_key, FUN = function(x) mean(x, na.rm = TRUE))
    d[[residual_col]] <- d$observed - d[[mean_col]]
  }

  d
}

residual_data <- add_split_mean_residuals(residual_data, split_var)

summarize_residuals <- function(d) {
  ok <- is.finite(d$observed) & is.finite(d$predicted)
  dd <- d[ok, , drop = FALSE]
  residual <- dd$observed - dd$predicted

  has_lm <- nrow(dd) >= 2 &&
    length(unique(dd$observed)) > 1 &&
    length(unique(dd$predicted)) > 1

  if (has_lm) {
    fit <- stats::lm(observed ~ predicted, data = dd)
    intercept <- unname(stats::coef(fit)[["(Intercept)"]])
    slope <- unname(stats::coef(fit)[["predicted"]])
    r2 <- summary(fit)$r.squared
    correlation <- stats::cor(dd$observed, dd$predicted)
  } else {
    intercept <- NA_real_
    slope <- NA_real_
    r2 <- NA_real_
    correlation <- NA_real_
  }

  data.frame(
    split_variable = split_var,
    split_value = as.character(dd[[split_var]][1]),
    model = dd$model[1],
    model_structure = dd$model_structure[1],
    hierarchy = dd$hierarchy[1],
    k_depth = dd$k_depth[1],
    response = dd$response[1],
    n = nrow(dd),
    n_stands = length(unique(dd$stand_id)),
    observed_mean = mean(dd$observed, na.rm = TRUE),
    predicted_mean = mean(dd$predicted, na.rm = TRUE),
    mean_residual = mean(residual, na.rm = TRUE),
    mean_abs_error = mean(abs(residual), na.rm = TRUE),
    rmse = sqrt(mean(residual^2, na.rm = TRUE)),
    residual_sd = stats::sd(residual, na.rm = TRUE),
    lm_intercept_observed_vs_predicted = intercept,
    lm_slope_observed_vs_predicted = slope,
    lm_r2_observed_vs_predicted = r2,
    correlation_observed_predicted = correlation,
    stringsAsFactors = FALSE
  )
}

split_groups <- split(
  residual_data,
  list(residual_data[[split_var]], residual_data$model, residual_data$response),
  drop = TRUE
)
residual_summary <- do.call(rbind, lapply(split_groups, summarize_residuals))
residual_summary <- residual_summary[order(
  residual_summary$split_value,
  residual_summary$response,
  residual_summary$model
), ]

split_counts <- aggregate(
  cbind(n_observations = source_row_id, n_stands = stand_id) ~ split_value,
  data = data.frame(
    split_value = residual_data[[split_var]],
    source_row_id = residual_data$source_row_id,
    stand_id = residual_data$stand_id
  ),
  FUN = function(x) length(unique(x))
)
split_counts <- split_counts[order(-split_counts$n_observations, split_counts$split_value), ]

write_output(residual_data, "01_residual_diagnostics_with_dominant_species.txt")
write_output(residual_summary, "02_residual_summary_by_dominant_species.txt")
write_output(split_counts, "00_dominant_species_split_counts.txt")

# ---- plotting-data ----
preferred_figure_responses <- c("befa.st", "beft.st")
figure_response_names <- preferred_figure_responses[
  preferred_figure_responses %in% unique(residual_data$response)
]

if (length(figure_response_names) == 0) {
  figure_response_names <- unique(residual_data$response)
}

response_labels <- response_label(figure_response_names)
response_labels[figure_response_names == "befa.st"] <- "BEFA"
response_labels[figure_response_names == "beft.st"] <- "BEFT"
response_plot_styles <- make_response_plot_styles(figure_response_names, labels = response_labels)
model_names <- model_info$model

make_mean_residuals_long <- function(d, mean_groups = c("h1", "h2")) {
  base_cols <- c(
    "model", "model_structure", "hierarchy", "k_depth", "family",
    "response", "row_id", "h1", "h2", "x", "observed", split_var
  )

  out <- lapply(mean_groups, function(mean_group) {
    mean_col <- paste0("mean_observed_", mean_group)
    residual_col <- paste0("residual_from_", mean_group, "_mean")

    if (!all(c(mean_group, mean_col, residual_col) %in% names(d))) {
      return(NULL)
    }

    data.frame(
      d[, base_cols, drop = FALSE],
      mean_group = mean_group,
      mean_level = d[[mean_group]],
      mean_observed = d[[mean_col]],
      residual_from_mean = d[[residual_col]],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

plot_split_residual_figures <- function(split_value) {
  d <- residual_data[residual_data[[split_var]] == split_value, , drop = FALSE]
  d <- d[d$response %in% figure_response_names, , drop = FALSE]

  if (nrow(d) < min_split_n) {
    return(character(0))
  }

  d$mean_observed_h1 <- d$mean_observed_h1_within_split
  d$residual_from_h1_mean <- d$residual_from_h1_mean_within_split
  d$mean_observed_h2 <- d$mean_observed_h2_within_split
  d$residual_from_h2_mean <- d$residual_from_h2_mean_within_split

  mean_residuals <- make_mean_residuals_long(d)
  split_stub <- file_stub(as.character(split_value))
  file_tag <- paste0(split_var, "_", split_stub)

  files <- character(0)

  file <- paste0("05_residual_vs_predicted_h1_h2_befa_beft_", file_tag, ".pdf")
  plot_grouped_residual_pdf(
    d,
    file.path(figure_dir, file),
    models = model_names,
    responses = figure_response_names,
    styles = response_plot_styles,
    x_var = "predicted",
    x_lab = "Predicted"
  )
  files <- c(files, file)

  file <- paste0("05_residual_vs_rsd_h1_h2_befa_beft_", file_tag, ".pdf")
  plot_grouped_residual_pdf(
    d,
    file.path(figure_dir, file),
    models = model_names,
    responses = figure_response_names,
    styles = response_plot_styles,
    x_var = "x",
    x_lab = "Relative stand density"
  )
  files <- c(files, file)

  mean_groups <- intersect(c("h1", "h2"), unique(mean_residuals$mean_group))

  for (mean_group in mean_groups) {
    file <- paste0("05_residual_vs_rsd_mean_", mean_group, "_befa_beft_", file_tag, ".pdf")
    plot_grouped_mean_residual_pdf(
      mean_residuals,
      file.path(figure_dir, file),
      models = model_names,
      responses = figure_response_names,
      styles = response_plot_styles,
      mean_group = mean_group
    )
    files <- c(files, file)
  }

  files
}

split_values <- split_counts$split_value[split_counts$n_observations >= min_split_n]
generated_figures <- unlist(lapply(split_values, plot_split_residual_figures), use.names = FALSE)

# ---- figure-annotation ----
figure_dictionary <- data.frame(
  file = generated_figures,
  source_data = "processed_data/bef_bayes_gamma/05_observed_vs_predicted.txt joined to processed_data/plot_biomass.txt",
  split_variable = split_var,
  purpose = "Residual diagnostic figure split by dominant species group, using within-split h1/h2 mean residuals.",
  stringsAsFactors = FALSE
)
write_tsv_table(figure_dictionary, file.path(figure_dir, "00_figure_dictionary.txt"))

writeLines(
  c(
    "BEF Bayesian gamma dominant-species residual figures",
    "",
    paste("Generated:", as.character(Sys.time())),
    paste("Source gamma output directory:", gamma_output_dir),
    paste("Source observation data:", source_data_file),
    paste("Analysis output directory:", analysis_dir),
    paste("Figure output directory:", figure_dir),
    paste("Split variable:", split_var),
    paste("Minimum observations per split:", min_split_n),
    "",
    "Notes:",
    "- The default split is dominant_genus, derived from species.1 in processed_data/plot_biomass.txt.",
    "- Set BEF_DOMINANT_SPLIT=dominant_species to split by the full dominant species name instead.",
    "- Set BEF_DOMINANT_SPLIT=tree_genus or tree_species to split by the sampled-tree taxon instead of the stand dominant species.",
    "- Residual-from-h1 and residual-from-h2 means are recomputed within each split group.",
    "- Numeric residual diagnostics and summaries are written to the analysis output directory."
  ),
  file.path(figure_dir, "README_FIGURES.txt")
)

writeLines(
  c(
    paste("Figure run generated at:", as.character(Sys.time())),
    paste("Project root:", project_root),
    paste("Split variable:", split_var),
    paste("Analysis output directory:", analysis_dir),
    paste("Figure output directory:", figure_dir),
    "",
    "Split counts:",
    paste(split_counts$split_value, split_counts$n_observations, sep = "\t"),
    "",
    "Figure files:",
    sort(generated_figures),
    "",
    "Figure metadata files:",
    c("00_figure_dictionary.txt", "README_FIGURES.txt", "FIGURE_RUN_SUMMARY.txt")
  ),
  file.path(figure_dir, "FIGURE_RUN_SUMMARY.txt")
)

figure_dir
