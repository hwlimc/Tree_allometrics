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

if (!dir.exists(file.path(project_root, "bayes_outputs"))) {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
}

setwd(project_root)

library(brms)
library(posterior)
library(loo)

script_args <- commandArgs(trailingOnly = TRUE)
max_step_arg <- grep("^--max-step=", script_args, value = TRUE)
max_step <- if (length(max_step_arg) > 0) {
  as.integer(sub("^--max-step=", "", max_step_arg[[1]]))
} else {
  Inf
}

if (is.na(max_step) || max_step < 1) {
  stop("--max-step must be a positive integer.")
}

run_id <- format(Sys.time(), "%Y%m%d%H%M%S")
run_generated_at <- as.character(Sys.time())
out_dir <- file.path("processed_data", "bef_bayes_dist_comparision")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- model-list ----
model_dir <- "bayes_outputs"
response_names <- c("befa.st", "befr.st")
family_pattern <- "gamma|lognormal|lnorm|stud|student|tdis|gaussian|normal|ndis"
model_file_pattern <- paste0("^xp_(.+)_(", family_pattern, ")_rsd_4-4k-99-15[.]rds$")

family_label <- function(family) {
  family <- tolower(family)
  ifelse(
    family %in% c("stud", "student", "tdis"),
    "student",
    ifelse(
      family %in% c("lognormal", "lnorm"),
      "lognormal",
      ifelse(family %in% c("gaussian", "normal", "ndis"), "gaussian", family)
    )
  )
}

model_structure_label <- function(structure_key) {
  out <- gsub("-", "_", structure_key)
  out <- gsub("sp_code", "sp", out)
  out[out == "none_k0"] <- "base_k0"
  out
}

discover_model_info <- function(model_dir) {
  files <- list.files(model_dir, pattern = "[.]rds$", full.names = TRUE)
  file_names <- basename(files)
  match_data <- regmatches(file_names, regexec(model_file_pattern, file_names))
  keep <- lengths(match_data) > 0

  if (!any(keep)) {
    stop("No Bayesian model .rds files found in ", model_dir)
  }

  match_data <- match_data[keep]
  structure_key <- vapply(match_data, `[[`, character(1), 2)
  family <- family_label(vapply(match_data, `[[`, character(1), 3))
  model_structure <- model_structure_label(structure_key)
  k_depth <- sub("^.*_(k[0-9]+)$", "\\1", structure_key)
  hierarchy <- sub("_[kK][0-9]+$", "", structure_key)
  hierarchy[hierarchy == "none"] <- "base"
  hierarchy <- gsub("-", "_", hierarchy)
  hierarchy <- gsub("sp_code", "sp", hierarchy)

  out <- data.frame(
    model = paste(model_structure, family, sep = "_"),
    model_structure = model_structure,
    hierarchy = hierarchy,
    k_depth = k_depth,
    family = family,
    model_file = files[keep],
    stringsAsFactors = FALSE
  )

  out[order(out$hierarchy, out$k_depth, out$model_structure, out$family), ]
}

model_info <- discover_model_info(model_dir)

missing_files <- model_info$model_file[!file.exists(model_info$model_file)]
if (length(missing_files) > 0) {
  stop("Missing model files:\n", paste(missing_files, collapse = "\n"))
}

model_info$model_time <- as.character(file.info(model_info$model_file)$mtime)
model_info$model_file_size <- file.info(model_info$model_file)$size
fits <- setNames(lapply(model_info$model_file, readRDS), model_info$model)

brms_family_label <- function(fit) {
  fam <- fit$family
  if (is.list(fam) && length(fam) > 0) {
    return(paste(vapply(fam, function(x) x$family, character(1)), collapse = ";"))
  }

  NA_character_
}

model_info$brms_family <- vapply(fits, brms_family_label, character(1))

# ---- helper-functions ----
write_txt <- function(x, file) {
  write.table(
    x,
    file.path(out_dir, file),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_run_summary <- function(last_step = NULL) {
  lines <- c(
    paste("Run ID:", run_id),
    paste("Run generated at:", run_generated_at),
    paste("Project root:", project_root),
    paste("Output directory:", out_dir),
    paste("Discovered model structures:", length(unique(model_info$model_structure))),
    paste("Discovered families:", paste(sort(unique(model_info$family)), collapse = ", "))
  )

  if (!is.null(last_step)) {
    lines <- c(lines, paste("Stopped after step:", last_step))
  }

  writeLines(
    c(lines, "", "Files:", list.files(out_dir)),
    file.path(out_dir, "RUN_SUMMARY.txt")
  )
}

finish_if_max_step <- function(step) {
  if (is.finite(max_step) && max_step <= step) {
    write_run_summary(last_step = step)
    message("Stopping after step ", step, " because --max-step=", max_step)
    quit(save = "no", status = 0)
  }
}

add_model <- function(df, model) {
  i <- match(model, model_info$model)
  data.frame(
    model = model,
    model_structure = model_info$model_structure[i],
    hierarchy = model_info$hierarchy[i],
    k_depth = model_info$k_depth[i],
    family = model_info$family[i],
    brms_family = model_info$brms_family[i],
    df,
    row.names = NULL
  )
}

response_info_for_fit <- function(fit) {
  data_names <- names(fit$data)

  if (all(c("y1m1", "y2") %in% data_names)) {
    return(data.frame(
      resp = c("y1m1", "y2"),
      response = response_names,
      transform = c("shift_y1m1", "identity"),
      jacobian_log_response = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    ))
  }

  if (all(c("z1", "z2") %in% data_names)) {
    return(data.frame(
      resp = c("z1", "z2"),
      response = response_names,
      transform = c("log_y1m1", "log_y2"),
      jacobian_log_response = c(TRUE, TRUE),
      stringsAsFactors = FALSE
    ))
  }

  stop("Unknown response columns in fit$data: ", paste(data_names, collapse = ", "))
}

resp_label <- function(resp_info) {
  resp_info$response[[1]]
}

backtransform <- function(resp_info, x) {
  switch(
    resp_info$transform[[1]],
    shift_y1m1 = x + 1,
    identity = x,
    log_y1m1 = exp(x) + 1,
    log_y2 = exp(x),
    stop("Unknown response transform: ", resp_info$transform[[1]])
  )
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_which_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_integer_ else which.max(replace(x, is.na(x), -Inf))
}

safe_which_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_integer_ else which.min(replace(x, is.na(x), Inf))
}

run_by_model <- function(fun) {
  do.call(rbind, Map(fun, fits, model_info$model))
}

run_by_model_response <- function(fun, model_index = seq_along(fits)) {
  if (length(model_index) == 0) {
    stop("No models selected.")
  }

  out <- list()
  k <- 1L

  for (i in model_index) {
    resp_info <- response_info_for_fit(fits[[i]])
    for (j in seq_len(nrow(resp_info))) {
      out[[k]] <- fun(fits[[i]], model_info$model[[i]], resp_info[j, , drop = FALSE])
      k <- k + 1L
    }
  }

  do.call(rbind, out)
}

chain_id_for_fit <- function(fit, ndraws) {
  nchains <- NULL

  if (inherits(fit$fit, "stanfit")) {
    nchains <- fit$fit@sim$chains
  } else if (is.function(fit$fit$num_chains)) {
    nchains <- fit$fit$num_chains()
  }

  if (is.null(nchains) || ndraws %% nchains != 0) {
    return(NULL)
  }

  rep(seq_len(nchains), each = ndraws / nchains)
}

compute_loo <- function(fit, resp_info, newdata) {
  resp <- resp_info$resp[[1]]
  log_lik_matrix <- log_lik(fit, resp = resp, newdata = newdata)

  if (anyNA(log_lik_matrix)) {
    stop("NA values found in log-likelihood for response ", resp)
  }

  if (isTRUE(resp_info$jacobian_log_response[[1]])) {
    log_lik_matrix <- sweep(log_lik_matrix, 2, newdata[[resp]], "-")
  }

  chain_id <- chain_id_for_fit(fit, nrow(log_lik_matrix))
  r_eff <- if (is.null(chain_id)) {
    NULL
  } else {
    relative_eff(exp(log_lik_matrix), chain_id = chain_id)
  }

  loo(log_lik_matrix, r_eff = r_eff)
}

# ---- manifest ----
manifest <- data.frame(
  run_id = run_id,
  model_info,
  n_rows = sapply(fits, function(x) nrow(x$data)),
  response_columns = vapply(fits, function(x) paste(names(x$data), collapse = ","), character(1)),
  stringsAsFactors = FALSE
)
write_txt(manifest, "00_model_manifest.txt")

# ---- step-1-mcmc ----
get_mcmc <- function(fit, model) {
  np <- nuts_params(fit)
  s <- summarise_draws(as_draws_df(fit))
  s <- s[!grepl("^Ymi_", s$variable), ]

  i_rhat <- safe_which_max(s$rhat)
  i_bulk <- safe_which_min(s$ess_bulk)
  i_tail <- safe_which_min(s$ess_tail)
  treedepth <- np$Value[np$Parameter == "treedepth__"]

  out <- data.frame(
    divergences = sum(np$Parameter == "divergent__" & np$Value == 1),
    max_treedepth = safe_max(treedepth),
    max_rhat = safe_max(s$rhat),
    max_rhat_parameter = if (is.na(i_rhat)) NA_character_ else s$variable[i_rhat],
    min_bulk_ess = safe_min(s$ess_bulk),
    min_bulk_ess_parameter = if (is.na(i_bulk)) NA_character_ else s$variable[i_bulk],
    min_tail_ess = safe_min(s$ess_tail),
    min_tail_ess_parameter = if (is.na(i_tail)) NA_character_ else s$variable[i_tail],
    stringsAsFactors = FALSE
  )

  add_model(out, model)
}

step1 <- run_by_model(get_mcmc)
step1$passes_mcmc_screen <- step1$divergences == 0 &
  (is.na(step1$max_rhat) | step1$max_rhat <= 1.01)
step1 <- step1[order(step1$model_structure, step1$family), ]
write_txt(step1, "01_mcmc_diagnostics.txt")
finish_if_max_step(1)

# ---- step-2-posterior-predictive-check ----
get_ppc_yrep <- function(fit, model, resp_info, ndraws = 1000) {
  resp <- resp_info$resp[[1]]
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  obs <- backtransform(resp_info, dat[[resp]])

  yrep <- posterior_predict(
    fit,
    newdata = dat,
    resp = resp,
    ndraws = ndraws
  )
  yrep <- backtransform(resp_info, yrep)

  yrep_median <- apply(yrep, 2, median)
  yrep_q05 <- apply(yrep, 2, quantile, 0.05)
  yrep_q95 <- apply(yrep, 2, quantile, 0.95)
  yrep_sd <- apply(yrep, 2, sd)

  pp_p_upper <- colMeans(sweep(yrep, 2, obs, `>=`))
  pp_p_two_tail <- 2 * pmin(pp_p_upper, 1 - pp_p_upper)

  out <- data.frame(
    response = resp_label(resp_info),
    row_id = seq_len(nrow(dat)),
    x = dat$x,
    observed = obs,
    yrep_median = yrep_median,
    yrep_q05 = yrep_q05,
    yrep_q95 = yrep_q95,
    residual_obs_minus_yrep_median = obs - yrep_median,
    std_residual = (obs - yrep_median) / yrep_sd,
    pp_p_upper = pp_p_upper,
    pp_p_two_tail = pp_p_two_tail,
    outside_yrep90 = obs < yrep_q05 | obs > yrep_q95,
    stringsAsFactors = FALSE
  )

  add_model(out, model)
}

summarize_ppcheck <- function(d) {
  data.frame(
    model = d$model[1],
    model_structure = d$model_structure[1],
    hierarchy = d$hierarchy[1],
    k_depth = d$k_depth[1],
    family = d$family[1],
    brms_family = d$brms_family[1],
    response = d$response[1],
    n = nrow(d),
    observed_mean = mean(d$observed, na.rm = TRUE),
    yrep_median_mean = mean(d$yrep_median, na.rm = TRUE),
    mean_residual = mean(d$residual_obs_minus_yrep_median, na.rm = TRUE),
    rmse = sqrt(mean(d$residual_obs_minus_yrep_median^2, na.rm = TRUE)),
    mean_abs_std_residual = mean(abs(d$std_residual), na.rm = TRUE),
    prop_outside_yrep90 = mean(d$outside_yrep90, na.rm = TRUE),
    prop_pp_twotail_lt_0.10 = mean(d$pp_p_two_tail < 0.10, na.rm = TRUE),
    prop_pp_twotail_lt_0.05 = mean(d$pp_p_two_tail < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

step2_yrep <- run_by_model_response(get_ppc_yrep)
write_txt(step2_yrep, "02_ppcheck_observed_yrep_pointwise.txt")

step2_summary <- do.call(
  rbind,
  lapply(
    split(step2_yrep, list(step2_yrep$model, step2_yrep$response), drop = TRUE),
    summarize_ppcheck
  )
)
step2_summary <- step2_summary[order(step2_summary$model_structure, step2_summary$response, step2_summary$family), ]
write_txt(step2_summary, "02_ppcheck_observed_yrep_summary.txt")
finish_if_max_step(2)

# ---- step-3-loo-distribution-comparison ----
get_loo <- function(fit, model, resp_info) {
  resp <- resp_info$resp[[1]]
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  lo <- compute_loo(fit, resp_info = resp_info, newdata = dat)
  est <- lo$estimates
  pk <- pareto_k_values(lo)

  out <- data.frame(
    response = resp_label(resp_info),
    n = nrow(dat),
    elpd_loo = est["elpd_loo", "Estimate"],
    elpd_loo_se = est["elpd_loo", "SE"],
    p_loo = est["p_loo", "Estimate"],
    looic = est["looic", "Estimate"],
    looic_se = est["looic", "SE"],
    max_pareto_k = max(pk, na.rm = TRUE),
    n_pareto_k_gt_0.7 = sum(pk > 0.7, na.rm = TRUE),
    n_pareto_k_gt_1.0 = sum(pk > 1.0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  add_model(out, model)
}

step3_loo <- run_by_model_response(get_loo)
structure_response_group <- interaction(
  step3_loo$model_structure,
  step3_loo$response,
  drop = TRUE
)
step3_loo$elpd_diff_from_structure_response_best <- ave(
  step3_loo$elpd_loo,
  structure_response_group,
  FUN = function(x) x - max(x, na.rm = TRUE)
)
step3_loo$looic_diff_from_structure_response_best <- ave(
  step3_loo$looic,
  structure_response_group,
  FUN = function(x) x - min(x, na.rm = TRUE)
)
step3_loo$rank_within_structure_response <- ave(
  -step3_loo$elpd_loo,
  structure_response_group,
  FUN = function(x) rank(x, ties.method = "first")
)
step3_loo$pareto_warning <- step3_loo$n_pareto_k_gt_0.7 > 0
step3_loo <- step3_loo[order(
  step3_loo$model_structure,
  step3_loo$response,
  step3_loo$rank_within_structure_response
), ]
write_txt(step3_loo, "03_loo_metrics_by_response.txt")

step3_distribution_comparison <- aggregate(
  cbind(elpd_loo, p_loo, looic, n, n_pareto_k_gt_0.7, n_pareto_k_gt_1.0) ~
    model + model_structure + hierarchy + k_depth + family + brms_family,
  data = step3_loo,
  FUN = sum
)
step3_distribution_comparison$max_pareto_k <- tapply(
  step3_loo$max_pareto_k,
  step3_loo$model,
  max,
  na.rm = TRUE
)[step3_distribution_comparison$model]

structure_group <- step3_distribution_comparison$model_structure
step3_distribution_comparison$elpd_diff_from_structure_best <- ave(
  step3_distribution_comparison$elpd_loo,
  structure_group,
  FUN = function(x) x - max(x, na.rm = TRUE)
)
step3_distribution_comparison$looic_diff_from_structure_best <- ave(
  step3_distribution_comparison$looic,
  structure_group,
  FUN = function(x) x - min(x, na.rm = TRUE)
)
step3_distribution_comparison$rank_within_structure <- ave(
  -step3_distribution_comparison$elpd_loo,
  structure_group,
  FUN = function(x) rank(x, ties.method = "first")
)
step3_distribution_comparison$pareto_warning <- step3_distribution_comparison$n_pareto_k_gt_0.7 > 0
step3_distribution_comparison <- step3_distribution_comparison[order(
  step3_distribution_comparison$model_structure,
  step3_distribution_comparison$rank_within_structure
), ]
write_txt(step3_distribution_comparison, "03_loo_distribution_comparison_by_structure.txt")

# ---- output-annotation ----
output_dictionary <- data.frame(
  file = c(
    "00_model_manifest.txt",
    "00_output_dictionary.txt",
    "01_mcmc_diagnostics.txt",
    "02_ppcheck_observed_yrep_pointwise.txt",
    "02_ppcheck_observed_yrep_summary.txt",
    "03_loo_metrics_by_response.txt",
    "03_loo_distribution_comparison_by_structure.txt",
    "README_OUTPUTS.txt",
    "RUN_SUMMARY.txt"
  ),
  step = c(
    "setup",
    "setup",
    "1",
    "2",
    "2",
    "3",
    "3",
    "summary",
    "summary"
  ),
  row_level = c(
    "one row per discovered model file",
    "one row per output file",
    "one row per model",
    "one row per model response observation",
    "one row per model response",
    "one row per model response",
    "one row per model, summed over responses",
    "text",
    "text"
  ),
  purpose = c(
    "Input model inventory and file metadata.",
    "Column and file descriptions for this distribution-comparison run.",
    "MCMC diagnostics used to screen whether each fitted distribution sampled adequately.",
    "Pointwise posterior predictive checks on the original response scale.",
    "Aggregate posterior predictive check metrics by model and response.",
    "Response-level LOO metrics with rank and differences within each model structure and response.",
    "Main distribution comparison: LOO totals ranked within each model structure.",
    "Human-readable summary of the distribution-comparison outputs.",
    "Run time, project path, output path, and output file listing."
  ),
  stringsAsFactors = FALSE
)
write_txt(output_dictionary, "00_output_dictionary.txt")

writeLines(
  c(
    "BEF Bayesian distribution-comparison outputs",
    "",
    paste("Generated:", run_generated_at),
    paste("Run ID:", run_id),
    paste("Project root:", project_root),
    paste("Output directory:", out_dir),
    paste("Discovered model structures:", length(unique(model_info$model_structure))),
    paste("Discovered families:", paste(sort(unique(model_info$family)), collapse = ", ")),
    "",
    "Step 1: MCMC diagnostics",
    "- Use 01_mcmc_diagnostics.txt to screen divergences, max Rhat, and ESS.",
    "",
    "Step 2: posterior predictive checks",
    "- Use 02_ppcheck_observed_yrep_summary.txt to compare yrep behavior across families within a model structure.",
    "- Values are back-transformed to the original BEF response scale.",
    "",
    "Step 3: LOO distribution comparison",
    "- Use 03_loo_distribution_comparison_by_structure.txt as the main ranking table.",
    "- rank_within_structure = 1 is the best distribution assumption for that model structure by summed elpd_loo.",
    "- elpd_diff_from_structure_best and looic_diff_from_structure_best are computed within each model structure.",
    "- Student/log-response fits are Jacobian-adjusted before LOO comparison when their data use z1/z2.",
    "- pareto_warning marks models with one or more Pareto k values above 0.7."
  ),
  file.path(out_dir, "README_OUTPUTS.txt")
)

write_run_summary(last_step = 3)
out_dir
