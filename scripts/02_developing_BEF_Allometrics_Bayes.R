
# ---- setup ----
setwd("/Users/hyli0001/wrd/b/Dynamic_allometrics/")
project_root <- getwd()
library(brms)
library(posterior)
library(loo)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_generated_at <- as.character(Sys.time())
out_dir <- file.path("bayes_outputs", "bef_bayes", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- model-list ----
model_info <- data.frame(
  model = c(
    "base_k0",
    "ftp_k0",
    "ftp_k1",
    "sp_k0",
    "sp_k1",
    "ftp_sp_k0",
    "ftp_sp_k1",
    "ftp_sp_k2"
  ),
  model_file = c(
    "bayes_outputs/xp_none_k0_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp_k0_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp_k1_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_sp_code_k0_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_sp_code_k1_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp-sp_code_k0_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp-sp_code_k1_gamma_rsd_4-4k-99-15.rds",
    "bayes_outputs/xp_ftp-sp_code_k2_gamma_rsd_4-4k-99-15.rds"
  ),
  stringsAsFactors = FALSE
)

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
  data.frame(model = model, df, row.names = NULL)
}

resp_label <- function(resp) {
  if (identical(resp, "y1m1")) "befa.st" else "befr.st"
}

backtransform <- function(resp, x) {
  if (identical(resp, "y1m1")) x + 1 else x
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

run_by_model_response <- function(fun) {
  do.call(
    rbind,
    unlist(
      Map(function(fit, model) {
        list(
          fun(fit, model, "y1m1"),
          fun(fit, model, "y2")
        )
      }, fits, model_info$model),
      recursive = FALSE
    )
  )
}

plot_grid <- function(n) {
  ncol <- min(4, n)
  nrow <- ceiling(n / ncol)
  c(nrow = nrow, ncol = ncol)
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

compute_loo <- function(fit, resp, newdata) {
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
  run_generated_at = run_generated_at,
  model_info,
  n_rows = sapply(fits, function(x) nrow(x$data)),
  stringsAsFactors = FALSE
)

write_txt(manifest, "00_model_manifest.txt")
manifest

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
step1

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
head(step2)

# ---- step-3-ppcheck-data ----
get_ppc_yrep <- function(fit, model, resp, ndraws = 1000) {
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  obs <- backtransform(resp, dat[[resp]])

  yrep <- posterior_predict(
    fit,
    newdata = dat,
    resp = resp,
    ndraws = ndraws
  )
  yrep <- backtransform(resp, yrep)

  yrep_median <- apply(yrep, 2, median)
  yrep_q05 <- apply(yrep, 2, quantile, 0.05)
  yrep_q95 <- apply(yrep, 2, quantile, 0.95)
  yrep_sd <- apply(yrep, 2, sd)

  pp_p_upper <- colMeans(sweep(yrep, 2, obs, `>=`))
  pp_p_two_tail <- 2 * pmin(pp_p_upper, 1 - pp_p_upper)

  out <- data.frame(
    response = resp_label(resp),
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
head(step3_yrep)

# ---- step-3-ppcheck-summary ----
summarize_ppcheck <- function(d) {
  data.frame(
    model = d$model[1],
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
step3_summary

# ---- step-3-ppcheck-density-plot ----
plot_ppc_density_one <- function(fit, model, resp, ndraws = 50) {
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  obs <- backtransform(resp, dat[[resp]])

  yrep <- posterior_predict(
    fit,
    newdata = dat,
    resp = resp,
    ndraws = ndraws
  )
  yrep <- backtransform(resp, yrep)

  obs_den <- density(obs, na.rm = TRUE)
  yrep_den <- lapply(seq_len(nrow(yrep)), function(i) density(yrep[i, ], na.rm = TRUE))
  y_lim <- range(c(obs_den$y, unlist(lapply(yrep_den, `[[`, "y"))), na.rm = TRUE)

  plot(
    obs_den,
    lwd = 2,
    col = "black",
    ylim = y_lim,
    xlab = resp_label(resp),
    main = paste(resp_label(resp), model, sep = " | ")
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

pdf(file.path(out_dir, "03_ppcheck_density_observed_yrep.pdf"), width = 12, height = 7)
old_par <- par(no.readonly = TRUE)
layout_dim <- plot_grid(length(fits))

for (resp in c("y1m1", "y2")) {
  par(
    mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (mod in names(fits)) {
    plot_ppc_density_one(fits[[mod]], mod, resp)
  }

  mtext(paste("Density PPC:", resp_label(resp)), outer = TRUE, cex = 1.2)
}

par(old_par)
dev.off()

# ---- step-4-observed-predicted-data ----
get_obs_pred <- function(fit, model, resp) {
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  pred <- fitted(fit, newdata = dat, resp = resp)[, "Estimate"]

  out <- data.frame(
    response = resp_label(resp),
    row_id = seq_len(nrow(dat)),
    x = dat$x,
    observed = backtransform(resp, dat[[resp]]),
    predicted = backtransform(resp, pred),
    stringsAsFactors = FALSE
  )

  out$residual <- out$observed - out$predicted
  add_model(out, model)
}

step4 <- run_by_model_response(get_obs_pred)
write_txt(step4, "04_observed_vs_predicted.txt")
head(step4)

# ---- step-4-observed-predicted-plot ----
pdf(file.path(out_dir, "04_observed_vs_predicted.pdf"), width = 12, height = 7)
old_par <- par(no.readonly = TRUE)
layout_dim <- plot_grid(length(fits))

for (resp in c("befa.st", "befr.st")) {
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
get_loo <- function(fit, model, resp) {
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  lo <- compute_loo(fit, resp = resp, newdata = dat)
  est <- lo$estimates
  pk <- pareto_k_values(lo)

  out <- data.frame(
    response = resp_label(resp),
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
step5

# ---- step-5-total-ranking ----
step5_total <- aggregate(
  cbind(elpd_loo, p_loo, looic, n, n_pareto_k_gt_0.7, n_pareto_k_gt_1.0) ~ model,
  data = step5,
  FUN = sum
)

step5_total$max_pareto_k <- tapply(step5$max_pareto_k, step5$model, max, na.rm = TRUE)[step5_total$model]
step5_total$elpd_diff_from_best <- step5_total$elpd_loo - max(step5_total$elpd_loo)
step5_total$looic_diff_from_best <- step5_total$looic - min(step5_total$looic)
step5_total <- step5_total[order(-step5_total$elpd_loo), ]

write_txt(step5_total, "05_loo_total_ranking.txt")
step5_total

# ---- step-5-response-rsd-data ----
get_response_curve <- function(fit, model, resp, n_grid = 100) {
  dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
  x_grid <- seq(min(dat$x, na.rm = TRUE), max(dat$x, na.rm = TRUE), length.out = n_grid)
  newdat <- data.frame(x = x_grid)

  if ("h1" %in% names(fit$data)) {
    newdat$h1 <- fit$data$h1[1]
  }
  if ("h2" %in% names(fit$data)) {
    newdat$h2 <- fit$data$h2[1]
  }

  pred <- fitted(
    fit,
    newdata = newdat,
    resp = resp,
    re_formula = NA,
    probs = c(0.05, 0.95)
  )

  obs <- add_model(data.frame(
    response = resp_label(resp),
    row_id = seq_len(nrow(dat)),
    x = dat$x,
    observed = backtransform(resp, dat[[resp]]),
    stringsAsFactors = FALSE
  ), model)

  curve <- add_model(data.frame(
    response = resp_label(resp),
    x = x_grid,
    predicted = backtransform(resp, pred[, "Estimate"]),
    pred_q05 = backtransform(resp, pred[, "Q5"]),
    pred_q95 = backtransform(resp, pred[, "Q95"]),
    stringsAsFactors = FALSE
  ), model)

  list(obs = obs, curve = curve)
}

curve_list <- unlist(
  Map(function(fit, model) {
    list(
      get_response_curve(fit, model, "y1m1"),
      get_response_curve(fit, model, "y2")
    )
  }, fits, model_info$model),
  recursive = FALSE
)

step5_obs <- do.call(rbind, lapply(curve_list, `[[`, "obs"))
step5_curve <- do.call(rbind, lapply(curve_list, `[[`, "curve"))

write_txt(step5_obs, "05_response_vs_rsd_observed.txt")
write_txt(step5_curve, "05_response_vs_rsd_predicted.txt")
head(step5_curve)

# ---- step-5-response-rsd-plot ----
pdf(file.path(out_dir, "05_response_vs_rsd.pdf"), width = 12, height = 7)
old_par <- par(no.readonly = TRUE)
layout_dim <- plot_grid(length(fits))

for (resp in c("befa.st", "befr.st")) {
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
