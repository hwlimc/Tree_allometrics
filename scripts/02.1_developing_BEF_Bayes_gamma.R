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
out_dir <- file.path("processed_data", "bef_bayes_gamma", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- gamma-model-list ----
model_dir <- "bayes_outputs"
target_family <- "gamma"
response_names <- c("befa.st", "befr.st")
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

plot_grid <- function(n) {
  ncol <- min(4, n)
  nrow <- ceiling(n / ncol)
  c(nrow = nrow, ncol = ncol)
}

plot_pdf_height <- function(n) {
  layout_dim <- plot_grid(n)
  max(7, 2.4 * layout_dim["nrow"])
}

plot_range <- function(x) {
  out <- range(x, na.rm = TRUE)

  if (!all(is.finite(out))) {
    return(c(0, 1))
  }

  if (diff(out) == 0) {
    pad <- if (out[1] == 0) 0.5 else abs(out[1]) * 0.05
    return(out + c(-pad, pad))
  }

  pad <- diff(out) * 0.04
  out + c(-pad, pad)
}

plot_group_label <- function(x, missing_label) {
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

# ---- manifest ----
manifest <- data.frame(
  run_id = run_id,
  model_info,
  n_rows = sapply(fits, function(x) nrow(x$data)),
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
write_txt(step1, "01_mcmc_diagnostics.txt")

# ---- step-2-posterior ----
get_posteriors <- function(fit, model) {
  s <- summarise_draws(as_draws_df(fit))
  s <- s[grepl("^(b_|sd_|shape_)", s$variable), ]

  keep <- intersect(
    c("variable", "mean", "median", "sd", "mad", "q5", "q95"),
    names(s)
  )

  add_model(s[, keep, drop = FALSE], model)
}

step2 <- run_by_model(get_posteriors)
write_txt(step2, "02_posterior_parameter_summary.txt")

# ---- step-3-ppcheck-data ----
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

step3_yrep <- run_by_model_response(get_ppc_yrep)
write_txt(step3_yrep, "03_ppcheck_observed_yrep_pointwise.txt")

# ---- step-3-ppcheck-summary ----
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

step3_summary <- do.call(
  rbind,
  lapply(split(step3_yrep, list(step3_yrep$model, step3_yrep$response), drop = TRUE), summarize_ppcheck)
)

step3_summary <- step3_summary[order(step3_summary$response, step3_summary$model), ]
write_txt(step3_summary, "03_ppcheck_observed_yrep_summary.txt")

# ---- step-3-ppcheck-density-plot ----
plot_ppc_density_one <- function(fit, model, resp_info, ndraws = 50) {
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

  obs_den <- density(obs, na.rm = TRUE)
  yrep_den <- lapply(seq_len(nrow(yrep)), function(i) density(yrep[i, ], na.rm = TRUE))
  y_lim <- range(c(obs_den$y, unlist(lapply(yrep_den, `[[`, "y"))), na.rm = TRUE)

  plot(
    obs_den,
    lwd = 2,
    col = "black",
    ylim = y_lim,
    xlab = resp_label(resp_info),
    main = paste(resp_label(resp_info), model, sep = " | ")
  )

  for (i in seq_along(yrep_den)) {
    lines(yrep_den[[i]], col = rgb(1, 0, 0, 0.15))
  }

  lines(obs_den, lwd = 2, col = "black")
  legend(
    "topright",
    legend = c("Observed", "yrep"),
    col = c("black", "red"),
    lwd = 2,
    bty = "n"
  )
}

layout_dim <- plot_grid(length(fits))
pdf(
  file.path(out_dir, "03_ppcheck_density_observed_yrep_gamma.pdf"),
  width = 12,
  height = plot_pdf_height(length(fits))
)
old_par <- par(no.readonly = TRUE)

for (response_name in response_names) {
  par(
    mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (i in seq_along(fits)) {
    resp_info <- response_info_for_fit(fits[[i]])
    resp_info <- resp_info[resp_info$response == response_name, , drop = FALSE]
    plot_ppc_density_one(fits[[i]], model_info$model[[i]], resp_info)
  }

  mtext(paste("Density PPC:", response_name, "| gamma"), outer = TRUE, cex = 1.2)
}

par(old_par)
dev.off()

# ---- step-4-observed-predicted-data ----
get_obs_pred <- function(fit, model, resp_info) {
  resp <- resp_info$resp[[1]]
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  pred <- fitted(fit, newdata = dat, resp = resp)[, "Estimate"]

  out <- data.frame(
    response = resp_label(resp_info),
    row_id = seq_len(nrow(dat)),
    x = dat$x,
    observed = backtransform(resp_info, dat[[resp]]),
    predicted = backtransform(resp_info, pred),
    stringsAsFactors = FALSE
  )

  out$residual <- out$observed - out$predicted
  add_model(out, model)
}

step4 <- run_by_model_response(get_obs_pred)
write_txt(step4, "04_observed_vs_predicted.txt")

# ---- step-4-observed-predicted-plot ----
layout_dim <- plot_grid(length(fits))
pdf(
  file.path(out_dir, "04_observed_vs_predicted.pdf"),
  width = 12,
  height = plot_pdf_height(length(fits))
)
old_par <- par(no.readonly = TRUE)

for (resp in response_names) {
  par(
    mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (mod in names(fits)) {
    d <- step4[step4$response == resp & step4$model == mod, ]
    lim <- range(c(d$observed, d$predicted), na.rm = TRUE)

    plot(
      d$observed,
      d$predicted,
      xlim = lim,
      ylim = lim,
      pch = 16,
      col = rgb(0, 0, 0, 0.35),
      xlab = "Observed",
      ylab = "Predicted",
      main = paste(resp, mod, sep = " | ")
    )

    abline(0, 1, col = "red", lwd = 2)
  }

  mtext(paste("Observed vs predicted:", resp), outer = TRUE, cex = 1.2)
}

par(old_par)
dev.off()

# ---- step-5-loo ----
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

step5 <- run_by_model_response(get_loo)
step5 <- step5[order(step5$response, -step5$elpd_loo), ]
write_txt(step5, "05_loo_metrics_by_response.txt")

# ---- step-5-total-ranking ----
step5_total <- aggregate(
  cbind(elpd_loo, p_loo, looic, n, n_pareto_k_gt_0.7, n_pareto_k_gt_1.0) ~
    model + model_structure + hierarchy + k_depth + family,
  data = step5,
  FUN = sum
)

step5_total$max_pareto_k <- tapply(step5$max_pareto_k, step5$model, max, na.rm = TRUE)[step5_total$model]
step5_total$elpd_diff_from_best <- step5_total$elpd_loo - max(step5_total$elpd_loo)
step5_total$looic_diff_from_best <- step5_total$looic - min(step5_total$looic)
step5_total <- step5_total[order(-step5_total$elpd_loo), ]
step5_total$rank <- seq_len(nrow(step5_total))

write_txt(step5_total, "05_loo_total_ranking.txt")

# ---- step-5-response-rsd-data ----
get_response_curve <- function(fit, model, resp_info, n_grid = 100) {
  resp <- resp_info$resp[[1]]
  grid <- prediction_grid_for_fit(fit, resp_info, n_grid = n_grid)
  dat <- grid$dat
  x_grid <- grid$x_grid

  pred <- fitted(
    fit,
    newdata = grid$newdata,
    resp = resp,
    re_formula = NA,
    probs = c(0.05, 0.95)
  )

  obs <- add_model(data.frame(
    response = resp_label(resp_info),
    row_id = seq_len(nrow(dat)),
    x = dat$x,
    observed = backtransform(resp_info, dat[[resp]]),
    stringsAsFactors = FALSE
  ), model)

  curve <- add_model(data.frame(
    response = resp_label(resp_info),
    x = x_grid,
    predicted = backtransform(resp_info, pred[, "Estimate"]),
    pred_q05 = backtransform(resp_info, pred[, "Q5"]),
    pred_q95 = backtransform(resp_info, pred[, "Q95"]),
    stringsAsFactors = FALSE
  ), model)

  list(obs = obs, curve = curve)
}

curve_list <- unlist(
  Map(function(fit, model) {
    resp_info <- response_info_for_fit(fit)
    lapply(seq_len(nrow(resp_info)), function(i) {
      get_response_curve(fit, model, resp_info[i, , drop = FALSE])
    })
  }, fits, model_info$model),
  recursive = FALSE
)

step5_obs <- do.call(rbind, lapply(curve_list, `[[`, "obs"))
step5_curve <- do.call(rbind, lapply(curve_list, `[[`, "curve"))

write_txt(step5_obs, "05_response_vs_rsd_observed.txt")
write_txt(step5_curve, "05_response_vs_rsd_predicted.txt")

# ---- step-5-response-rsd-h1-h2-data ----
get_grouped_response_curve <- function(fit, model, resp_info, n_grid = 100) {
  resp <- resp_info$resp[[1]]
  dat <- fit$data[!is.na(fit$data[[resp]]) & !is.na(fit$data$x), , drop = FALSE]
  group_vars <- intersect(c("h1", "h2"), names(dat))

  obs <- add_model(data.frame(
    response = resp_label(resp_info),
    row_id = seq_len(nrow(dat)),
    h1 = if ("h1" %in% names(dat)) plot_group_label(dat$h1, "missing h1") else rep("all", nrow(dat)),
    h2 = if ("h2" %in% names(dat)) plot_group_label(dat$h2, "missing h2") else rep("all", nrow(dat)),
    x = dat$x,
    observed = backtransform(resp_info, dat[[resp]]),
    stringsAsFactors = FALSE
  ), model)

  if (length(group_vars) == 0) {
    combos <- data.frame(.all = 1L)
    combo_rows <- list(seq_len(nrow(dat)))
  } else {
    complete <- complete.cases(dat[, group_vars, drop = FALSE])
    combos <- unique(dat[complete, group_vars, drop = FALSE])

    if (nrow(combos) == 0) {
      stop("No complete h1/h2 groups available for ", model, " ", resp_label(resp_info))
    }

    combos <- combos[do.call(order, lapply(combos, function(x) as.character(x))), , drop = FALSE]
    combo_rows <- lapply(seq_len(nrow(combos)), function(i) {
      keep <- complete

      for (v in group_vars) {
        keep <- keep & dat[[v]] == combos[[v]][i]
      }

      which(keep)
    })
  }

  curve_parts <- lapply(seq_along(combo_rows), function(i) {
    rows <- combo_rows[[i]]
    dat_i <- dat[rows, , drop = FALSE]
    x_min <- min(dat_i$x, na.rm = TRUE)
    x_max <- max(dat_i$x, na.rm = TRUE)
    x_grid <- if (isTRUE(all.equal(x_min, x_max))) x_min else seq(x_min, x_max, length.out = n_grid)
    newdat <- data.frame(x = x_grid)

    for (v in group_vars) {
      newdat[[v]] <- combos[[v]][rep(i, length(x_grid))]
    }

    pred <- fitted(
      fit,
      newdata = newdat,
      resp = resp,
      re_formula = NULL,
      probs = c(0.05, 0.95)
    )

    data.frame(
      response = resp_label(resp_info),
      h1 = if ("h1" %in% group_vars) plot_group_label(combos$h1[i], "missing h1") else "all",
      h2 = if ("h2" %in% group_vars) plot_group_label(combos$h2[i], "missing h2") else "all",
      x_min_observed = x_min,
      x_max_observed = x_max,
      x = x_grid,
      predicted = backtransform(resp_info, pred[, "Estimate"]),
      pred_q05 = backtransform(resp_info, pred[, "Q5"]),
      pred_q95 = backtransform(resp_info, pred[, "Q95"]),
      stringsAsFactors = FALSE
    )
  })

  list(
    obs = obs,
    curve = add_model(do.call(rbind, curve_parts), model)
  )
}

grouped_curve_list <- unlist(
  Map(function(fit, model) {
    resp_info <- response_info_for_fit(fit)
    lapply(seq_len(nrow(resp_info)), function(i) {
      get_grouped_response_curve(fit, model, resp_info[i, , drop = FALSE])
    })
  }, fits, model_info$model),
  recursive = FALSE
)

step5_grouped_obs <- do.call(rbind, lapply(grouped_curve_list, `[[`, "obs"))
step5_grouped_curve <- do.call(rbind, lapply(grouped_curve_list, `[[`, "curve"))

write_txt(step5_grouped_obs, "05_response_vs_rsd_h1_h2_observed.txt")
write_txt(step5_grouped_curve, "05_response_vs_rsd_h1_h2_predicted.txt")

# ---- step-6-publication-figure-data ----
# These files are small, annotated derivatives for publication figures. They keep
# posterior uncertainty needed for plotting without writing full draw matrices.
best_gamma_model <- step5_total$model[[1]]

add_best_gamma_flag <- function(df) {
  df$is_best_gamma_model <- df$model == best_gamma_model
  df
}

get_publication_curve_summary <- function(fit, model, resp_info) {
  resp <- resp_info$resp[[1]]
  grid <- prediction_grid_for_fit(fit, resp_info)

  epred <- posterior_epred(
    fit,
    newdata = grid$newdata,
    resp = resp,
    re_formula = NA
  )
  epred <- backtransform(resp_info, epred)

  q <- t(apply(epred, 2, quantile, probs = publication_curve_probs, na.rm = TRUE))
  colnames(q) <- paste0("pred_", publication_curve_prob_names)

  add_model(data.frame(
    response = resp_label(resp_info),
    x = grid$x_grid,
    predicted_mean = colMeans(epred, na.rm = TRUE),
    predicted_median = apply(epred, 2, median, na.rm = TRUE),
    predicted_sd = apply(epred, 2, sd, na.rm = TRUE),
    q,
    stringsAsFactors = FALSE
  ), model)
}

get_publication_curve_draws <- function(fit, model, resp_info) {
  resp <- resp_info$resp[[1]]
  grid <- prediction_grid_for_fit(fit, resp_info)

  epred <- posterior_epred(
    fit,
    newdata = grid$newdata,
    resp = resp,
    re_formula = NA,
    ndraws = publication_curve_draws_n
  )
  epred <- backtransform(resp_info, epred)

  add_model(data.frame(
    response = resp_label(resp_info),
    draw_id = rep(seq_len(nrow(epred)), each = length(grid$x_grid)),
    x = rep(grid$x_grid, times = nrow(epred)),
    predicted = as.vector(t(epred)),
    stringsAsFactors = FALSE
  ), model)
}

get_publication_parameter_summary <- function(fit, model) {
  draws <- as_draws_df(fit)
  variables <- names(draws)
  variables <- variables[grepl("^(b_|sd_|shape_)", variables)]

  out <- do.call(rbind, lapply(variables, function(variable) {
    x <- draws[[variable]]
    data.frame(
      variable = variable,
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      q025 = unname(quantile(x, 0.025, na.rm = TRUE)),
      q05 = unname(quantile(x, 0.05, na.rm = TRUE)),
      q95 = unname(quantile(x, 0.95, na.rm = TRUE)),
      q975 = unname(quantile(x, 0.975, na.rm = TRUE)),
      prob_gt_0 = mean(x > 0, na.rm = TRUE),
      prob_lt_0 = mean(x < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  add_model(out, model)
}

step6_curve_summary <- do.call(
  rbind,
  lapply(seq_along(fits), function(i) {
    resp_info <- response_info_for_fit(fits[[i]])
    do.call(rbind, lapply(seq_len(nrow(resp_info)), function(j) {
      get_publication_curve_summary(fits[[i]], model_info$model[[i]], resp_info[j, , drop = FALSE])
    }))
  })
)
step6_curve_summary <- add_best_gamma_flag(step6_curve_summary)
write_txt(step6_curve_summary, "06_publication_response_curve_summary_gamma.txt")

step6_curve_draws <- do.call(
  rbind,
  lapply(seq_along(fits), function(i) {
    resp_info <- response_info_for_fit(fits[[i]])
    do.call(rbind, lapply(seq_len(nrow(resp_info)), function(j) {
      get_publication_curve_draws(fits[[i]], model_info$model[[i]], resp_info[j, , drop = FALSE])
    }))
  })
)
step6_curve_draws <- add_best_gamma_flag(step6_curve_draws)
write_txt(step6_curve_draws, "06_publication_response_curve_draws_gamma.txt")

step6_parameter_summary <- do.call(
  rbind,
  lapply(seq_along(fits), function(i) {
    get_publication_parameter_summary(fits[[i]], model_info$model[[i]])
  })
)
step6_parameter_summary <- add_best_gamma_flag(step6_parameter_summary)
write_txt(step6_parameter_summary, "06_publication_parameter_summary_gamma.txt")

# ---- step-5-response-rsd-plot ----
layout_dim <- plot_grid(length(fits))
pdf(
  file.path(out_dir, "05_response_vs_rsd.pdf"),
  width = 12,
  height = plot_pdf_height(length(fits))
)
old_par <- par(no.readonly = TRUE)

for (resp in response_names) {
  par(
    mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (mod in names(fits)) {
    obs_d <- step5_obs[step5_obs$response == resp & step5_obs$model == mod, ]
    cur_d <- step5_curve[step5_curve$response == resp & step5_curve$model == mod, ]
    ord <- order(cur_d$x)
    ylim <- range(c(obs_d$observed, cur_d$pred_q05, cur_d$pred_q95), na.rm = TRUE)

    plot(
      obs_d$x,
      obs_d$observed,
      pch = 16,
      col = rgb(0, 0, 0, 0.35),
      xlab = "Relative stand density (rsd)",
      ylab = resp,
      main = paste(resp, mod, sep = " | "),
      ylim = ylim
    )

    lines(cur_d$x[ord], cur_d$predicted[ord], col = "red", lwd = 2)
    lines(cur_d$x[ord], cur_d$pred_q05[ord], col = "red", lty = 2)
    lines(cur_d$x[ord], cur_d$pred_q95[ord], col = "red", lty = 2)
  }

  mtext(paste("Response vs rsd:", resp), outer = TRUE, cex = 1.2)
}

par(old_par)
dev.off()

# ---- step-5-response-rsd-h1-h2-plot ----
max_h1_panels <- max(vapply(
  split(step5_grouped_obs, list(step5_grouped_obs$response, step5_grouped_obs$model), drop = TRUE),
  function(d) length(unique(d$h1)),
  integer(1)
))

pdf(
  file.path(out_dir, "05_response_vs_rsd_h1_h2.pdf"),
  width = 12,
  height = plot_pdf_height(max_h1_panels)
)
old_par <- par(no.readonly = TRUE)

for (resp in response_names) {
  for (mod in names(fits)) {
    obs_d <- step5_grouped_obs[step5_grouped_obs$response == resp & step5_grouped_obs$model == mod, ]
    cur_d <- step5_grouped_curve[step5_grouped_curve$response == resp & step5_grouped_curve$model == mod, ]
    h1_levels <- sort(unique(obs_d$h1))
    h2_levels <- sort(unique(obs_d$h2))
    h2_cols <- setNames(hcl.colors(length(h2_levels), palette = "Dark 3"), h2_levels)
    layout_dim <- plot_grid(length(h1_levels))
    layout_slots <- prod(layout_dim)
    xlim <- plot_range(c(obs_d$x, cur_d$x))
    ylim <- plot_range(c(obs_d$observed, cur_d$pred_q05, cur_d$pred_q95))

    par(
      mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
      mar = c(4, 4, 3, 1),
      oma = c(0, 0, 2, 0)
    )

    for (panel_i in seq_along(h1_levels)) {
      h1_level <- h1_levels[[panel_i]]

      plot(
        NA,
        NA,
        xlim = xlim,
        ylim = ylim,
        xlab = "Relative stand density (rsd)",
        ylab = resp,
        main = paste("h1:", h1_level)
      )

      for (h2_level in h2_levels) {
        obs_i <- obs_d[obs_d$h1 == h1_level & obs_d$h2 == h2_level, ]
        cur_i <- cur_d[cur_d$h1 == h1_level & cur_d$h2 == h2_level, ]
        col_i <- h2_cols[[h2_level]]

        if (nrow(obs_i) > 0) {
          points(obs_i$x, obs_i$observed, pch = 16, col = adjustcolor(col_i, alpha.f = 0.45))
        }

        if (nrow(cur_i) > 0) {
          ord <- order(cur_i$x)
          lines(cur_i$x[ord], cur_i$predicted[ord], col = col_i, lwd = 2)
          lines(cur_i$x[ord], cur_i$pred_q05[ord], col = adjustcolor(col_i, alpha.f = 0.70), lty = 2)
          lines(cur_i$x[ord], cur_i$pred_q95[ord], col = adjustcolor(col_i, alpha.f = 0.70), lty = 2)
        }
      }

      if (panel_i == 1 && any(h2_levels != "all")) {
        legend_ncol <- if (length(h2_levels) > 16) 3 else if (length(h2_levels) > 8) 2 else 1
        legend_cex <- if (length(h2_levels) > 16) 0.45 else if (length(h2_levels) > 8) 0.55 else 0.70

        legend(
          "topright",
          legend = h2_levels,
          col = h2_cols,
          pch = 16,
          lty = 1,
          lwd = 2,
          cex = legend_cex,
          ncol = legend_ncol,
          bty = "n",
          title = "h2"
        )
      }
    }

    if (layout_slots > length(h1_levels)) {
      for (i in seq_len(layout_slots - length(h1_levels))) {
        plot.new()
      }
    }

    mtext(paste("Response vs rsd by h1/h2:", resp, "|", mod), outer = TRUE, cex = 1.2)
  }
}

par(old_par)
dev.off()

# ---- output-annotation ----
output_dictionary <- data.frame(
  file = c(
    "00_model_manifest.txt",
    "01_mcmc_diagnostics.txt",
    "02_posterior_parameter_summary.txt",
    "03_ppcheck_observed_yrep_pointwise.txt",
    "03_ppcheck_observed_yrep_summary.txt",
    "03_ppcheck_density_observed_yrep_gamma.pdf",
    "04_observed_vs_predicted.txt",
    "04_observed_vs_predicted.pdf",
    "05_loo_metrics_by_response.txt",
    "05_loo_total_ranking.txt",
    "05_response_vs_rsd_observed.txt",
    "05_response_vs_rsd_predicted.txt",
    "05_response_vs_rsd.pdf",
    "05_response_vs_rsd_h1_h2_observed.txt",
    "05_response_vs_rsd_h1_h2_predicted.txt",
    "05_response_vs_rsd_h1_h2.pdf",
    "06_publication_response_curve_summary_gamma.txt",
    "06_publication_response_curve_draws_gamma.txt",
    "06_publication_parameter_summary_gamma.txt",
    "RUN_SUMMARY.txt"
  ),
  models = c(
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "gamma models only",
    "run metadata"
  ),
  row_level = c(
    "one row per model",
    "one row per model",
    "one row per model parameter",
    "one row per model response observation",
    "one row per model response",
    "plot",
    "one row per model response observation",
    "plot",
    "one row per gamma model response",
    "one row per gamma model",
    "one row per model response observation",
    "one row per model response x-grid point",
    "plot",
    "one row per model response observation",
    "one row per model response h1/h2 x-grid point",
    "plot",
    "one row per gamma model response x-grid point",
    "one row per sampled posterior draw and x-grid point",
    "one row per gamma model parameter",
    "text"
  ),
  response_scale = c(
    "metadata",
    "diagnostic scale",
    "parameter scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "LOO on gamma response scale",
    "summed LOO on gamma response scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "original response scale",
    "parameter scale",
    "metadata"
  ),
  purpose = c(
    "Input gamma model inventory and file metadata.",
    "MCMC diagnostics for screening divergent or weakly mixed gamma fits.",
    "Posterior summaries for population, group SD, and gamma shape parameters.",
    "Pointwise posterior predictive checks using observed values and yrep summaries.",
    "Aggregate posterior predictive check metrics for model screening.",
    "Density posterior predictive checks for gamma models.",
    "Observed versus posterior fitted medians for calibration plots.",
    "Quick base-R observed versus predicted diagnostic plot.",
    "Response-level LOO metrics for comparing gamma hierarchical structures.",
    "Total LOO ranking for comparing gamma hierarchical structures.",
    "Observed points for response versus relative stand density figures.",
    "Posterior fitted response curves with 90 percent intervals.",
    "Quick base-R response versus relative stand density plot.",
    "Observed points with h1 panel labels and h2 color labels.",
    "Posterior fitted response curves by h1/h2, clipped to each observed h1/h2 rsd range.",
    "Response versus relative stand density panels by h1 with h2-colored points and curves.",
    "Publication-ready gamma response curves with multiple credible intervals.",
    "Sampled gamma posterior expected curves for spaghetti or ribbon diagnostics.",
    "Publication-ready gamma parameter summaries with posterior probabilities.",
    "Run time, project path, output path, and output file listing."
  ),
  stringsAsFactors = FALSE
)
write_txt(output_dictionary, "00_output_dictionary.txt")

writeLines(
  c(
    "BEF Bayesian gamma-only post-processing outputs",
    "",
    paste("Generated:", run_generated_at),
    paste("Run ID:", run_id),
    paste("Project root:", project_root),
    paste("Output directory:", out_dir),
    paste("Discovered gamma models:", paste(model_info$model, collapse = ", ")),
    paste("Best gamma model:", best_gamma_model),
    "",
    "Scale notes:",
    "- Only gamma model files are loaded by this workflow.",
    "- Figure datasets are on the original response scale: befa.st and befr.st.",
    "- LOO compares hierarchical structures within the gamma distribution.",
    "- Publication curve files use posterior expected responses, not full posterior predictive yrep draws.",
    "- 06_publication_response_curve_draws_gamma.txt stores a compact draw sample for plotting, not all posterior draws.",
    "- 05_response_vs_rsd_h1_h2.pdf uses conditional h1/h2 predictions where those grouping levels exist.",
    "- Grouped h1/h2 prediction lines are clipped to the observed rsd range for each h1/h2 combination.",
    "",
    "Recommended publication inputs:",
    "- Main fitted-curve figure: 06_publication_response_curve_summary_gamma.txt.",
    "- Optional spaghetti curves: 06_publication_response_curve_draws_gamma.txt.",
    "- Grouped h1/h2 diagnostic figure: 05_response_vs_rsd_h1_h2.pdf.",
    "- Gamma hierarchy ranking: 05_loo_total_ranking.txt.",
    "- Diagnostics and screening: 01_mcmc_diagnostics.txt and 03_ppcheck_observed_yrep_summary.txt.",
    "- Column and file descriptions: 00_output_dictionary.txt."
  ),
  file.path(out_dir, "README_OUTPUTS.txt")
)

# ---- run-summary ----
writeLines(
  c(
    paste("Run ID:", run_id),
    paste("Run generated at:", run_generated_at),
    paste("Project root:", project_root),
    paste("Output directory:", out_dir),
    "",
    "Files:",
    list.files(out_dir)
  ),
  file.path(out_dir, "RUN_SUMMARY.txt")
)

out_dir
