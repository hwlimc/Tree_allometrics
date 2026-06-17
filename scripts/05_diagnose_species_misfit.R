## Diagnose species-level misfit for saved brms BEF models.
##
## Example:
## Rscript --vanilla scripts/05_diagnose_species_misfit.R \
##   bayes_outputs/exp_decay_log_partial_pool_ftp_sp_code_kdepth0_rsd_student_4chn_4000itr_4cor_0.99del_15depth.rds \
##   sp_code
##
## Optional arguments:
##   1 fit_file
##   2 group_col, usually sp_code; if absent, h2/h1/grp are tried
##   3 out_dir
##   4 ndraws

suppressPackageStartupMessages({
  library(brms)
})

if (basename(getwd()) == "scripts") {
  setwd("..")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    "Usage: Rscript --vanilla scripts/05_diagnose_species_misfit.R ",
    "<fit_file> [group_col] [out_dir] [ndraws]"
  )
}

fit_file <- args[1]
group_col_arg <- if (length(args) >= 2) args[2] else "sp_code"
out_dir <- if (length(args) >= 3) {
  args[3]
} else {
  file.path(
    dirname(fit_file),
    paste0(tools::file_path_sans_ext(basename(fit_file)), "_misfit")
  )
}
ndraws <- if (length(args) >= 4) as.integer(args[4]) else 1000L

if (!file.exists(fit_file)) {
  stop("Fit file does not exist: ", fit_file)
}
if (is.na(ndraws) || ndraws < 1) {
  stop("ndraws must be a positive integer.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(fit_file)
dat <- fit$data

group_candidates <- unique(c(group_col_arg, "sp_code", "h2", "h1", "grp"))
group_col <- group_candidates[group_candidates %in% names(dat)][1]
if (is.na(group_col)) {
  stop(
    "Could not find a grouping column. Tried: ",
    paste(group_candidates, collapse = ", ")
  )
}

if (!"x" %in% names(dat)) {
  stop("Expected fit$data to contain column 'x'.")
}

response_vars <- intersect(c("y1", "y2"), names(dat))
if (length(response_vars) == 0) {
  stop("Expected fit$data to contain y1 and/or y2.")
}

log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(mean(exp(x - max_x)))
}

as_draw_matrix <- function(x) {
  x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  x
}

diagnose_response <- function(response) {
  obs_id <- which(!is.na(dat[[response]]))
  if (length(obs_id) == 0) {
    return(NULL)
  }

  newdata <- dat[obs_id, , drop = FALSE]
  observed <- newdata[[response]]
  group <- as.character(newdata[[group_col]])

  message("Predicting ", response, " for ", length(obs_id), " observed rows")
  epred <- as_draw_matrix(
    posterior_epred(
      fit,
      newdata = newdata,
      resp = response,
      re_formula = NULL,
      ndraws = ndraws
    )
  )
  yrep <- as_draw_matrix(
    posterior_predict(
      fit,
      newdata = newdata,
      resp = response,
      re_formula = NULL,
      ndraws = ndraws
    )
  )

  pred_median <- apply(epred, 2, stats::median)
  pred_q05 <- apply(epred, 2, stats::quantile, probs = 0.05)
  pred_q95 <- apply(epred, 2, stats::quantile, probs = 0.95)
  yrep_q05 <- apply(yrep, 2, stats::quantile, probs = 0.05)
  yrep_q95 <- apply(yrep, 2, stats::quantile, probs = 0.95)
  yrep_sd <- apply(yrep, 2, stats::sd)

  p_upper <- colMeans(sweep(yrep, 2, observed, `>=`))
  p_two_tail <- 2 * pmin(p_upper, 1 - p_upper)
  residual <- observed - pred_median

  pointwise <- data.frame(
    response = response,
    row_id = obs_id,
    group = group,
    x = newdata$x,
    observed = observed,
    epred_median = pred_median,
    epred_q05 = pred_q05,
    epred_q95 = pred_q95,
    pred_q05 = yrep_q05,
    pred_q95 = yrep_q95,
    residual = residual,
    std_residual = residual / yrep_sd,
    pp_p_upper = p_upper,
    pp_p_two_tail = p_two_tail,
    outside_pred90 = observed < yrep_q05 | observed > yrep_q95,
    stringsAsFactors = FALSE
  )

  group_levels <- sort(unique(group))
  summary_rows <- lapply(group_levels, function(g) {
    i <- which(group == g)
    bias_draws <- rowMeans(sweep(epred[, i, drop = FALSE], 2, observed[i], `-`))
    bias_draws <- -bias_draws

    data.frame(
      response = response,
      group = g,
      n = length(i),
      x_min = min(newdata$x[i], na.rm = TRUE),
      x_max = max(newdata$x[i], na.rm = TRUE),
      x_range = diff(range(newdata$x[i], na.rm = TRUE)),
      obs_mean = mean(observed[i], na.rm = TRUE),
      pred_mean = mean(pred_median[i], na.rm = TRUE),
      mean_residual = mean(residual[i], na.rm = TRUE),
      rmse = sqrt(mean(residual[i]^2, na.rm = TRUE)),
      mean_abs_std_residual = mean(abs(pointwise$std_residual[i]), na.rm = TRUE),
      prop_outside_pred90 = mean(pointwise$outside_pred90[i], na.rm = TRUE),
      prop_pp_twotail_lt_0.10 = mean(pointwise$pp_p_two_tail[i] < 0.10, na.rm = TRUE),
      prop_pp_twotail_lt_0.05 = mean(pointwise$pp_p_two_tail[i] < 0.05, na.rm = TRUE),
      bias_q05 = stats::quantile(bias_draws, 0.05, na.rm = TRUE),
      bias_q50 = stats::quantile(bias_draws, 0.50, na.rm = TRUE),
      bias_q95 = stats::quantile(bias_draws, 0.95, na.rm = TRUE),
      prob_bias_positive = mean(bias_draws > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary_rows)

  summary$systematic_bias_90 <- summary$bias_q05 > 0 | summary$bias_q95 < 0
  summary$many_outside_pred90 <- summary$n >= 5 & summary$prop_outside_pred90 > 0.20
  summary$many_extreme_pp <- summary$n >= 5 & summary$prop_pp_twotail_lt_0.10 > 0.20
  summary$misfit_flag <- with(
    summary,
    systematic_bias_90 | many_outside_pred90 | many_extreme_pp
  )

  list(pointwise = pointwise, summary = summary)
}

results <- lapply(response_vars, diagnose_response)
names(results) <- response_vars
results <- Filter(Negate(is.null), results)

pointwise <- do.call(rbind, lapply(results, `[[`, "pointwise"))
summary <- do.call(rbind, lapply(results, `[[`, "summary"))
summary <- summary[order(summary$misfit_flag, summary$response, -summary$rmse), ]
summary <- summary[order(!summary$misfit_flag, summary$response, -summary$rmse), ]

write.csv(pointwise, file.path(out_dir, "species_misfit_pointwise.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "species_misfit_summary.csv"), row.names = FALSE)
saveRDS(
  list(
    fit_file = fit_file,
    group_col = group_col,
    ndraws = ndraws,
    pointwise = pointwise,
    summary = summary
  ),
  file.path(out_dir, "species_misfit.rds")
)

pdf(file.path(out_dir, "species_misfit_diagnostics.pdf"), width = 10, height = 7)
for (response in unique(summary$response)) {
  s <- summary[summary$response == response, , drop = FALSE]
  s <- s[order(-s$rmse), ]
  top <- head(s, 20)
  old_mar <- par("mar")
  par(mar = c(8, 5, 4, 2))
  barplot(
    top$rmse,
    names.arg = top$group,
    las = 2,
    ylab = "RMSE",
    main = paste(response, "species ranked by RMSE")
  )
  par(mar = old_mar)

  flagged <- s[s$misfit_flag, , drop = FALSE]
  if (nrow(flagged) > 0) {
    flagged <- head(flagged[order(abs(flagged$bias_q50), decreasing = TRUE), ], 20)
    plot(
      seq_len(nrow(flagged)),
      flagged$bias_q50,
      ylim = range(c(flagged$bias_q05, flagged$bias_q95), na.rm = TRUE),
      xaxt = "n",
      xlab = "",
      ylab = "Mean observed - expected",
      main = paste(response, "flagged species mean residual")
    )
    arrows(
      seq_len(nrow(flagged)), flagged$bias_q05,
      seq_len(nrow(flagged)), flagged$bias_q95,
      angle = 90,
      code = 3,
      length = 0.04
    )
    abline(h = 0, lty = 2)
    axis(1, at = seq_len(nrow(flagged)), labels = flagged$group, las = 2)
  }
}
dev.off()

message("\nSpecies-level misfit summary:")
print(summary[summary$misfit_flag, ], row.names = FALSE)

message("\nWrote outputs to: ", normalizePath(out_dir))
