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

run_id <- format(Sys.time(), "%Y%m%d")
run_generated_at <- as.character(Sys.time())
out_dir <- file.path("processed_data", "bef_bayes_gamma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- gamma-model-list ----
model_dir <- "bayes_outputs"
target_family <- "gamma"
response_names <- c("befa.st", "befr.st")
total_response_name <- "beft.st"
publication_curve_grid_n <- 100
publication_curve_draws_n <- 200
publication_curve_probs <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
publication_curve_prob_names <- c("q025", "q05", "q10", "q25", "q50", "q75", "q90", "q95", "q975")
model_file_pattern <- "^xp_(.+)_gamma_rsd_4-4k-99-15[.]rds$"

model_structure_label <- function(structure_key) {
  out <- gsub("-", "_", structure_key)
  out <- gsub("sp_code", "sp", out)
  out[out == "none_k0"] <- "base_k0"
  out
}

discover_gamma_model_info <- function(model_dir) {
  files <- list.files(model_dir, pattern = "[.]rds$", full.names = TRUE)
  file_names <- basename(files)
  match_data <- regmatches(file_names, regexec(model_file_pattern, file_names))
  keep <- lengths(match_data) > 0

  if (!any(keep)) {
    stop("No gamma Bayesian model .rds files found in ", model_dir)
  }

  match_data <- match_data[keep]
  structure_key <- vapply(match_data, `[[`, character(1), 2)
  model_structure <- model_structure_label(structure_key)
  k_depth <- sub("^.*_(k[0-9]+)$", "\\1", structure_key)
  hierarchy <- sub("_[kK][0-9]+$", "", structure_key)
  hierarchy[hierarchy == "none"] <- "base"
  hierarchy <- gsub("-", "_", hierarchy)
  hierarchy <- gsub("sp_code", "sp", hierarchy)

  out <- data.frame(
    model = paste(model_structure, target_family, sep = "_"),
    model_structure = model_structure,
    hierarchy = hierarchy,
    k_depth = k_depth,
    family = target_family,
    model_file = files[keep],
    stringsAsFactors = FALSE
  )

  out[order(out$hierarchy, out$k_depth, out$model_structure), ]
}

model_info <- discover_gamma_model_info(model_dir)

missing_files <- model_info$model_file[!file.exists(model_info$model_file)]
if (length(missing_files) > 0) {
  stop("Missing model files:\n", paste(missing_files, collapse = "\n"))
}

model_info$model_time <- as.character(file.info(model_info$model_file)$mtime)
model_info$model_file_size <- file.info(model_info$model_file)$size

fits <- setNames(lapply(model_info$model_file, readRDS), model_info$model)

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

add_model <- function(df, model) {
  i <- match(model, model_info$model)
  data.frame(
    model = model,
    model_structure = model_info$model_structure[i],
    hierarchy = model_info$hierarchy[i],
    k_depth = model_info$k_depth[i],
    family = model_info$family[i],
    df,
    row.names = NULL
  )
}

response_info_for_fit <- function(fit) {
  data_names <- names(fit$data)

  if (!all(c("y1m1", "y2") %in% data_names)) {
    stop(
      "Gamma fits should contain response columns y1m1 and y2. Found: ",
      paste(data_names, collapse = ", ")
    )
  }

  data.frame(
    resp = c("y1m1", "y2"),
    response = response_names,
    transform = c("shift_y1m1", "identity"),
    stringsAsFactors = FALSE
  )
}

resp_label <- function(resp_info) {
  resp_info$response[[1]]
}

backtransform <- function(resp_info, x) {
  switch(
    resp_info$transform[[1]],
    shift_y1m1 = x + 1,
    identity = x,
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
    stop("No gamma models selected.")
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

group_label <- function(x, missing_label) {
  out <- as.character(x)
  out[is.na(out)] <- missing_label
  out
}

prediction_grid_for_fit <- function(fit, resp_info, n_grid = publication_curve_grid_n) {
  resp <- resp_info$resp[[1]]
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  x_grid <- seq(min(dat$x, na.rm = TRUE), max(dat$x, na.rm = TRUE), length.out = n_grid)
  newdat <- data.frame(x = x_grid)

  if ("h1" %in% names(fit$data)) {
    newdat$h1 <- fit$data$h1[1]
  }
  if ("h2" %in% names(fit$data)) {
    newdat$h2 <- fit$data$h2[1]
  }

  list(dat = dat, x_grid = x_grid, newdata = newdat)
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

  chain_id <- chain_id_for_fit(fit, nrow(log_lik_matrix))
  r_eff <- if (is.null(chain_id)) {
    NULL
  } else {
    relative_eff(exp(log_lik_matrix), chain_id = chain_id)
  }

  loo(log_lik_matrix, r_eff = r_eff)
}


##############################################################################
##############################################################################
## ---------------------------------------------------------------------------
##        | Posterior for   | Main object          | Main question
## ---------------------------------------------------------------------------
## Step 1 | Check MCMC      | Rhat, ESS, divergence| Did sampling work?
## Step 2 | Simulate data   | yrep vs y            | Can model mimic data?
## Step 3 | Compare models  | LOO, Pareto k        | Which model predicts best?
## Step 4 | Summarize params| b_, sd_, shape_      | What was estimated?
## Step 5 | Fit observations| predicted vs observed| Where does fit miss?
## Step 6 | Make plot data  | curves, residuals    | What to plot/report?
## ---------------------------------------------------------------------------
##############################################################################
##############################################################################


# ---- model information ----
model_info <- data.frame(
  run_id = run_id,
  model_info,
  n_rows = sapply(fits, function(x) nrow(x$data)),
  stringsAsFactors = FALSE)
write_txt(model_info, "00_model_info.txt")

# ---- Step 1: mcmc check ----
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
write_txt(step1, "01_mcmc_diagnostics.txt")



# ---- Step 2: posteriror predictive check ----
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

step2_yrep <- run_by_model_response(get_ppc_yrep)
write_txt(step2_yrep, "02_ppcheck_obs_vs_yrep.txt")

summarize_ppcheck <- function(d) {
  data.frame(
    model = d$model[1],
    model_structure = d$model_structure[1],
    hierarchy = d$hierarchy[1],
    k_depth = d$k_depth[1],
    family = d$family[1],
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

step2_summary <- do.call(
  rbind,
  lapply(split(step2_yrep, list(step2_yrep$model, step2_yrep$response), drop = TRUE), summarize_ppcheck)
)

step2_summary <- step2_summary[order(step2_summary$response, step2_summary$model), ]
write_txt(step2_summary, "02_ppcheck_summary.txt")


# ---- step-3-loo-model-comparison ----
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
step3_loo <- step3_loo[order(step3_loo$response, -step3_loo$elpd_loo), ]
write_txt(step3_loo, "03_loo_metrics_by_response.txt")

# ---- step-3-total-ranking ----
step3_total <- aggregate(
  cbind(elpd_loo, p_loo, looic, n, n_pareto_k_gt_0.7, n_pareto_k_gt_1.0) ~
    model + model_structure + hierarchy + k_depth + family,
  data = step3_loo,
  FUN = sum
)

step3_total$max_pareto_k <- tapply(step3_loo$max_pareto_k, step3_loo$model, max, na.rm = TRUE)[step3_total$model]
step3_total$elpd_diff_from_best <- step3_total$elpd_loo - max(step3_total$elpd_loo)
step3_total$looic_diff_from_best <- step3_total$looic - min(step3_total$looic)
step3_total <- step3_total[order(-step3_total$elpd_loo), ]
step3_total$rank <- seq_len(nrow(step3_total))

write_txt(step3_total, "03_loo_total_ranking.txt")


# ---- step-4-posterior-parameter-summary ----
get_posteriors <- function(fit, model) {
  s <- summarise_draws(as_draws_df(fit))
  s <- s[grepl("^(b_|sd_|shape_)", s$variable), ]

  keep <- intersect(
    c("variable", "mean", "median", "sd", "mad", "q5", "q95"),
    names(s)
  )

  add_model(s[, keep, drop = FALSE], model)
}

step4_parameters <- run_by_model(get_posteriors)
write_txt(step4_parameters, "04_posterior_parameter_summary.txt")


# ---- Step 5: Observed vs predicted ----
get_observed_fit <- function(fit, model, resp_info) {
  resp <- resp_info$resp[[1]]
  dat <- fit$data[!is.na(fit$data[[resp]]) & !is.na(fit$data$x), , drop = FALSE]
  pred <- fitted(
    fit,
    newdata = dat,
    resp = resp,
    re_formula = NULL
  )[, "Estimate"]

  observed <- backtransform(resp_info, dat[[resp]])
  predicted <- backtransform(resp_info, pred)

  add_model(data.frame(
    response = resp_label(resp_info),
    row_id = seq_len(nrow(dat)),
    h1 = if ("h1" %in% names(dat)) group_label(dat$h1, "missing h1") else rep("all", nrow(dat)),
    h2 = if ("h2" %in% names(dat)) group_label(dat$h2, "missing h2") else rep("all", nrow(dat)),
    x = dat$x,
    observed = observed,
    predicted = predicted,
    residual = observed - predicted,
    stringsAsFactors = FALSE
  ), model)
}

get_observed_total_fit <- function(fit, model) {
  resp_info <- response_info_for_fit(fit)
  befa_info <- resp_info[resp_info$response == "befa.st", , drop = FALSE]
  befr_info <- resp_info[resp_info$response == "befr.st", , drop = FALSE]
  dat <- fit$data[!is.na(fit$data$y1m1) & !is.na(fit$data$y2) & !is.na(fit$data$x), , drop = FALSE]

  pred_befa <- fitted(
    fit,
    newdata = dat,
    resp = befa_info$resp[[1]],
    re_formula = NULL
  )[, "Estimate"]
  pred_befr <- fitted(
    fit,
    newdata = dat,
    resp = befr_info$resp[[1]],
    re_formula = NULL
  )[, "Estimate"]

  observed <- backtransform(befa_info, dat$y1m1) + backtransform(befr_info, dat$y2)
  predicted <- backtransform(befa_info, pred_befa) + backtransform(befr_info, pred_befr)

  add_model(data.frame(
    response = total_response_name,
    row_id = seq_len(nrow(dat)),
    h1 = if ("h1" %in% names(dat)) group_label(dat$h1, "missing h1") else rep("all", nrow(dat)),
    h2 = if ("h2" %in% names(dat)) group_label(dat$h2, "missing h2") else rep("all", nrow(dat)),
    x = dat$x,
    observed = observed,
    predicted = predicted,
    residual = observed - predicted,
    stringsAsFactors = FALSE
  ), model)
}

add_mean_residual_columns <- function(d, mean_groups = c("h1", "h2")) {
  for (mean_group in mean_groups) {
    if (!mean_group %in% names(d)) {
      next
    }

    group_key <- interaction(
      d$model,
      d$response,
      d[[mean_group]],
      drop = TRUE
    )
    mean_col <- paste0("mean_observed_", mean_group)
    residual_col <- paste0("residual_from_", mean_group, "_mean")

    d[[mean_col]] <- ave(d$observed, group_key, FUN = function(x) mean(x, na.rm = TRUE))
    d[[residual_col]] <- d$observed - d[[mean_col]]
  }

  d
}

summarize_observed_fit <- function(d) {
  ok <- is.finite(d$observed) & is.finite(d$predicted)
  dd <- d[ok, , drop = FALSE]
  residual <- dd$observed - dd$predicted

  has_lm <- nrow(dd) >= 2 &&
    length(unique(dd$observed)) > 1 &&
    length(unique(dd$predicted)) > 1

  if (has_lm) {
    fit <- lm(observed ~ predicted, data = dd)
    intercept <- unname(coef(fit)[["(Intercept)"]])
    slope <- unname(coef(fit)[["predicted"]])
    r2 <- summary(fit)$r.squared
    correlation <- cor(dd$observed, dd$predicted)
  } else {
    intercept <- NA_real_
    slope <- NA_real_
    r2 <- NA_real_
    correlation <- NA_real_
  }

  data.frame(
    model = d$model[1],
    model_structure = d$model_structure[1],
    hierarchy = d$hierarchy[1],
    k_depth = d$k_depth[1],
    family = d$family[1],
    response = d$response[1],
    n = nrow(dd),
    observed_mean = mean(dd$observed, na.rm = TRUE),
    predicted_mean = mean(dd$predicted, na.rm = TRUE),
    mean_residual = mean(residual, na.rm = TRUE),
    mean_abs_error = mean(abs(residual), na.rm = TRUE),
    rmse = sqrt(mean(residual^2, na.rm = TRUE)),
    residual_sd = sd(residual, na.rm = TRUE),
    lm_intercept_observed_vs_predicted = intercept,
    lm_slope_observed_vs_predicted = slope,
    lm_r2_observed_vs_predicted = r2,
    correlation_observed_predicted = correlation,
    stringsAsFactors = FALSE
  )
}

observed_fit_list <- unlist(
  Map(function(fit, model) {
    resp_info <- response_info_for_fit(fit)
    c(lapply(seq_len(nrow(resp_info)), function(i) {
      get_observed_fit(fit, model, resp_info[i, , drop = FALSE])
    }), list(get_observed_total_fit(fit, model)))
  }, fits, model_info$model),
  recursive = FALSE
)

step5_observed_fit <- do.call(rbind, observed_fit_list)
step5_observed_fit <- add_mean_residual_columns(step5_observed_fit)
step5_observed_fit <- step5_observed_fit[, c(
  "model", "model_structure", "hierarchy", "k_depth", "family",
  "response", "row_id", "h1", "h2", "x",
  "observed", "predicted", "residual",
  "mean_observed_h1", "residual_from_h1_mean",
  "mean_observed_h2", "residual_from_h2_mean"
)]

step5_fit_summary <- do.call(
  rbind,
  lapply(
    split(step5_observed_fit, list(step5_observed_fit$model, step5_observed_fit$response), drop = TRUE),
    summarize_observed_fit
  )
)
step5_fit_summary <- step5_fit_summary[order(step5_fit_summary$response, step5_fit_summary$model), ]

write_txt(step5_observed_fit, "05_observed_vs_predicted.txt")
write_txt(step5_fit_summary, "05_observed_vs_predicted_summary.txt")

# # # ---- step-6-publication-figure-data ----
# # # These files are small, annotated derivatives for publication figures. They keep
# # # posterior uncertainty needed for plotting without writing full draw matrices.
# # best_gamma_model <- step3_total$model[[1]]

# # add_best_gamma_flag <- function(df) {
  # # df$is_best_gamma_model <- df$model == best_gamma_model
  # # df
# # }

# # get_publication_curve_summary <- function(fit, model, resp_info) {
  # # resp <- resp_info$resp[[1]]
  # # grid <- prediction_grid_for_fit(fit, resp_info)

  # # epred <- posterior_epred(
    # # fit,
    # # newdata = grid$newdata,
    # # resp = resp,
    # # re_formula = NA
  # # )
  # # epred <- backtransform(resp_info, epred)

  # # q <- t(apply(epred, 2, quantile, probs = publication_curve_probs, na.rm = TRUE))
  # # colnames(q) <- paste0("pred_", publication_curve_prob_names)

  # # add_model(data.frame(
    # # response = resp_label(resp_info),
    # # x = grid$x_grid,
    # # predicted_mean = colMeans(epred, na.rm = TRUE),
    # # predicted_median = apply(epred, 2, median, na.rm = TRUE),
    # # predicted_sd = apply(epred, 2, sd, na.rm = TRUE),
    # # q,
    # # stringsAsFactors = FALSE
  # # ), model)
# # }

# # get_publication_curve_draws <- function(fit, model, resp_info) {
  # # resp <- resp_info$resp[[1]]
  # # grid <- prediction_grid_for_fit(fit, resp_info)

  # # epred <- posterior_epred(
    # # fit,
    # # newdata = grid$newdata,
    # # resp = resp,
    # # re_formula = NA,
    # # ndraws = publication_curve_draws_n
  # # )
  # # epred <- backtransform(resp_info, epred)

  # # add_model(data.frame(
    # # response = resp_label(resp_info),
    # # draw_id = rep(seq_len(nrow(epred)), each = length(grid$x_grid)),
    # # x = rep(grid$x_grid, times = nrow(epred)),
    # # predicted = as.vector(t(epred)),
    # # stringsAsFactors = FALSE
  # # ), model)
# # }

# # get_publication_parameter_summary <- function(fit, model) {
  # # draws <- as_draws_df(fit)
  # # variables <- names(draws)
  # # variables <- variables[grepl("^(b_|sd_|shape_)", variables)]

  # # out <- do.call(rbind, lapply(variables, function(variable) {
    # # x <- draws[[variable]]
    # # data.frame(
      # # variable = variable,
      # # mean = mean(x, na.rm = TRUE),
      # # median = median(x, na.rm = TRUE),
      # # sd = sd(x, na.rm = TRUE),
      # # q025 = unname(quantile(x, 0.025, na.rm = TRUE)),
      # # q05 = unname(quantile(x, 0.05, na.rm = TRUE)),
      # # q95 = unname(quantile(x, 0.95, na.rm = TRUE)),
      # # q975 = unname(quantile(x, 0.975, na.rm = TRUE)),
      # # prob_gt_0 = mean(x > 0, na.rm = TRUE),
      # # prob_lt_0 = mean(x < 0, na.rm = TRUE),
      # # stringsAsFactors = FALSE
    # # )
  # # }))

  # # add_model(out, model)
# # }

# # step6_curve_summary <- do.call(
  # # rbind,
  # # lapply(seq_along(fits), function(i) {
    # # resp_info <- response_info_for_fit(fits[[i]])
    # # do.call(rbind, lapply(seq_len(nrow(resp_info)), function(j) {
      # # get_publication_curve_summary(fits[[i]], model_info$model[[i]], resp_info[j, , drop = FALSE])
    # # }))
  # # })
# # )
# # step6_curve_summary <- add_best_gamma_flag(step6_curve_summary)
# # write_txt(step6_curve_summary, "06_publication_response_curve_summary_gamma.txt")

# # step6_curve_draws <- do.call(
  # # rbind,
  # # lapply(seq_along(fits), function(i) {
    # # resp_info <- response_info_for_fit(fits[[i]])
    # # do.call(rbind, lapply(seq_len(nrow(resp_info)), function(j) {
      # # get_publication_curve_draws(fits[[i]], model_info$model[[i]], resp_info[j, , drop = FALSE])
    # # }))
  # # })
# # )
# # step6_curve_draws <- add_best_gamma_flag(step6_curve_draws)
# # write_txt(step6_curve_draws, "06_publication_response_curve_draws_gamma.txt")

# # step6_parameter_summary <- do.call(
  # # rbind,
  # # lapply(seq_along(fits), function(i) {
    # # get_publication_parameter_summary(fits[[i]], model_info$model[[i]])
  # # })
# # )
# # step6_parameter_summary <- add_best_gamma_flag(step6_parameter_summary)
# # write_txt(step6_parameter_summary, "06_publication_parameter_summary_gamma.txt")

# # # ---- output-annotation ----
# # output_dictionary <- data.frame(
  # # file = c(
    # # "00_model_manifest.txt",
    # # "00_output_dictionary.txt",
    # # "01_mcmc_diagnostics.txt",
    # # "02.1_ppcheck_obs_vs_yrep.txt",
    # # "02.2_ppcheck_summary.txt",
    # # "03.1_loo_metrics_by_response.txt",
    # # "03.2_loo_total_ranking.txt",
    # # "04_posterior_parameter_summary.txt",
    # # "05.1_observed_fit_diagnostics_h1_h2.txt",
    # # "05.2_observed_fit_summary.txt",
    # # "06_publication_response_curve_summary_gamma.txt",
    # # "06_publication_response_curve_draws_gamma.txt",
    # # "06_publication_parameter_summary_gamma.txt",
    # # "README_OUTPUTS.txt",
    # # "RUN_SUMMARY.txt"
  # # ),
  # # models = c(
    # # "gamma models only",
    # # "run metadata",
    # # rep("gamma models only", 11),
    # # "run metadata",
    # # "run metadata"
  # # ),
  # # row_level = c(
    # # "one row per model",
    # # "one row per output file",
    # # "one row per model",
    # # "one row per model response observation",
    # # "one row per model response",
    # # "one row per gamma model response",
    # # "one row per gamma model",
    # # "one row per model parameter",
    # # "one row per model response observation",
    # # "one row per model response",
    # # "one row per gamma model response x-grid point",
    # # "one row per sampled posterior draw and x-grid point",
    # # "one row per gamma model parameter",
    # # "text",
    # # "text"
  # # ),
  # # response_scale = c(
    # # "metadata",
    # # "metadata",
    # # "diagnostic scale",
    # # "original response scale",
    # # "original response scale",
    # # "LOO on gamma response scale",
    # # "summed LOO on gamma response scale",
    # # "parameter scale",
    # # "original response scale",
    # # "original response scale",
    # # "original response scale",
    # # "original response scale",
    # # "parameter scale",
    # # "metadata",
    # # "metadata"
  # # ),
  # # purpose = c(
    # # "Input gamma model inventory and file metadata.",
    # # "Column and file descriptions for numeric outputs from this script.",
    # # "MCMC diagnostics for screening divergent or weakly mixed gamma fits.",
    # # "Pointwise posterior predictive checks using observed values and yrep summaries.",
    # # "Aggregate posterior predictive check metrics for model screening.",
    # # "Response-level LOO metrics for comparing gamma hierarchical structures.",
    # # "Total LOO ranking for comparing gamma hierarchical structures.",
    # # "Posterior summaries for population, group SD, and gamma shape parameters.",
    # # "Observed, conditional predicted, residual, and h1/h2 mean-residual diagnostics in one long table.",
    # # "Observed-versus-predicted fit summary with RMSE, MAE, residual bias, linear-fit slope/intercept, R2, and correlation.",
    # # "Publication-ready gamma response curves with multiple credible intervals.",
    # # "Sampled gamma posterior expected curves for spaghetti or ribbon diagnostics.",
    # # "Publication-ready gamma parameter summaries with posterior probabilities.",
    # # "Human-readable summary of numeric outputs and figure-script handoff.",
    # # "Run time, project path, output path, and output file listing."
  # # ),
  # # stringsAsFactors = FALSE
# # )
# # write_txt(output_dictionary, "00_output_dictionary.txt")

# # writeLines(
  # # c(
    # # "BEF Bayesian gamma-only numeric outputs",
    # # "",
    # # paste("Generated:", run_generated_at),
    # # paste("Run ID:", run_id),
    # # paste("Project root:", project_root),
    # # paste("Output directory:", out_dir),
    # # paste("Discovered gamma models:", paste(model_info$model, collapse = ", ")),
    # # paste("Best gamma model:", best_gamma_model),
    # # "",
    # # "Scale notes:",
    # # "- This script writes tabular numeric outputs and metadata only; it does not write PDFs.",
    # # "- Figure-ready datasets are on the original response scale: befa.st, befr.st, and derived beft.st.",
    # # "- 05.1_observed_fit_diagnostics_h1_h2.txt keeps observed, predicted, residual, and h1/h2 mean-residual values in one long table.",
    # # "- 05.2_observed_fit_summary.txt summarizes observed-vs-predicted fit with RMSE, MAE, bias, linear-fit R2, slope, intercept, and correlation.",
    # # "- Figure scripts can overlay befa.st and beft.st in the same h1 panels; befr.st is retained in tables as the component used to compute beft.st.",
    # # "- LOO compares hierarchical structures within the gamma distribution.",
    # # "- Publication curve files use posterior expected responses, not full posterior predictive yrep draws.",
    # # "- 06_publication_response_curve_draws_gamma.txt stores a compact draw sample for plotting, not all posterior draws.",
    # # "- Mean-baseline residual columns use observed group means, not posterior predictions: residual = observed - mean(observed).",
    # # "",
    # # "Figure generation:",
    # # "- Render PDFs with: Rscript scripts/03_plot_BEF_Bayes_gamma_figures.R <output directory>",
    # # "- If <output directory> is omitted, the figure script uses the latest processed_data/bef_bayes_gamma run.",
    # # "- The PPC density figure is generated by the figure script from model .rds files listed in 00_model_manifest.txt.",
    # # "- Figure file descriptions are written by the figure script to 00_figure_dictionary.txt.",
    # # "",
    # # "Recommended publication inputs:",
    # # "- Main fitted-curve data: 06_publication_response_curve_summary_gamma.txt.",
    # # "- Optional spaghetti curve data: 06_publication_response_curve_draws_gamma.txt.",
    # # "- Observed fit and residual diagnostics: 05.1_observed_fit_diagnostics_h1_h2.txt.",
    # # "- Observed fit summary: 05.2_observed_fit_summary.txt.",
    # # "- Gamma hierarchy ranking: 03.2_loo_total_ranking.txt.",
    # # "- Parameter summaries: 04_posterior_parameter_summary.txt.",
    # # "- Diagnostics and screening: 01_mcmc_diagnostics.txt and 02.2_ppcheck_summary.txt.",
    # # "- Column and file descriptions: 00_output_dictionary.txt."
  # # ),
  # # file.path(out_dir, "README_OUTPUTS.txt")
# # )

# # # ---- run-summary ----
# # writeLines(
  # # c(
    # # paste("Run ID:", run_id),
    # # paste("Run generated at:", run_generated_at),
    # # paste("Project root:", project_root),
    # # paste("Output directory:", out_dir),
    # # "",
    # # "Numeric files:",
    # # output_dictionary$file
  # # ),
  # # file.path(out_dir, "RUN_SUMMARY.txt")
# # )

# # out_dir
