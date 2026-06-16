## Compare one-level and two-level Bayesian allometry models.
##
## Primary criterion: PSIS-LOO expected log predictive density (ELPD).
## Because y2 was fitted with mi(), compare y1 on all rows and y2 on rows
## where y2 is observed.

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

if (basename(getwd()) == "scripts") {
  setwd("..")
}

model_files <- c(
  ftp_sp = "bayes_outputs/exp_decay_ftp_sp_code_rsd_4chn_4000itr_4cor_0.99del_15depth.rds",
  ftp    = "bayes_outputs/exp_decay_ftp_rsd_4chn_4000itr_4cor_0.99del_15depth.rds",
  sp_c   = "bayes_outputs/exp_decay_sp_code_rsd_4chn_4000itr_4cor_0.99del_15depth.rds",
  pft_sp = "bayes_outputs/exp_decay_PFT_sp_code_rsd_4chn_4000itr_4cor_0.99del_15depth.rds",
  pft    = "bayes_outputs/exp_decay_PFT_rsd_4chn_4000itr_4cor_0.99del_15depth.rds"
)

model_info <- data.frame(
  model = names(model_files),
  hierarchy = c("ftp + sp_code", "ftp", "sp_code", "PFT + sp_code", "PFT"),
  hierarchy_depth = c(2L, 1L, 1L, 2L, 1L),
  stringsAsFactors = FALSE
)

out_dir <- "bayes_outputs/model_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

force_recompute <- FALSE

loo_cache_file <- function(model, response) {
  file.path(out_dir, paste0("loo_", model, "_", response, ".rds"))
}

read_cached_loo <- function(model, response) {
  cache_file <- loo_cache_file(model, response)
  if (force_recompute || !file.exists(cache_file)) {
    return(NULL)
  }

  cache <- readRDS(cache_file)
  model_file <- model_files[[model]]
  if (!identical(cache$model_file, model_file)) {
    return(NULL)
  }
  if (!identical(cache$model_mtime, file.info(model_file)$mtime)) {
    return(NULL)
  }
  if (!identical(cache$method, "log_lik_matrix")) {
    return(NULL)
  }
  cache$loo
}

write_cached_loo <- function(model, response, loo_object) {
  model_file <- model_files[[model]]
  saveRDS(
    list(
      model = model,
      response = response,
      model_file = model_file,
      model_mtime = file.info(model_file)$mtime,
      method = "log_lik_matrix",
      loo = loo_object
    ),
    loo_cache_file(model, response)
  )
}

chain_id_for_fit <- function(fit, ndraws) {
  if (inherits(fit$fit, "stanfit")) {
    nchains <- fit$fit@sim$chains
  } else if (is.function(fit$fit$num_chains)) {
    nchains <- fit$fit$num_chains()
  } else {
    nchains <- 1L
    warning("Could not determine chain count; assuming independent draws.")
  }

  if (ndraws %% nchains != 0) {
    warning("Number of draws is not divisible by chains; assuming independent draws.")
    return(NULL)
  }

  rep(seq_len(nchains), each = ndraws / nchains)
}

compute_loo <- function(fit, model, response, newdata = NULL) {
  cached <- read_cached_loo(model, response)
  if (!is.null(cached)) {
    message("Using cached LOO: ", model, " / ", response)
    return(cached)
  }

  message("Computing LOO: ", model, " / ", response)
  args <- list(object = fit, resp = response)
  if (!is.null(newdata)) {
    args$newdata <- newdata
  }

  log_lik_matrix <- do.call(brms::log_lik, args)
  if (anyNA(log_lik_matrix)) {
    stop("NA values found in log-likelihood for ", model, " / ", response)
  }
  chain_id <- chain_id_for_fit(fit, nrow(log_lik_matrix))
  r_eff <- if (is.null(chain_id)) {
    NULL
  } else {
    loo::relative_eff(exp(log_lik_matrix), chain_id = chain_id)
  }
  loo_object <- loo::loo(log_lik_matrix, r_eff = r_eff)

  write_cached_loo(model, response, loo_object)
  loo_object
}

loo_metrics <- function(loo_object, model, response) {
  estimates <- loo_object$estimates
  pareto_k <- loo::pareto_k_values(loo_object)
  data.frame(
    model = model,
    response = response,
    n_units = nrow(loo_object$pointwise),
    elpd_loo = estimates["elpd_loo", "Estimate"],
    elpd_loo_se = estimates["elpd_loo", "SE"],
    p_loo = estimates["p_loo", "Estimate"],
    looic = estimates["looic", "Estimate"],
    looic_se = estimates["looic", "SE"],
    max_pareto_k = max(pareto_k, na.rm = TRUE),
    n_pareto_k_gt_0.7 = sum(pareto_k > 0.7, na.rm = TRUE),
    n_pareto_k_gt_1 = sum(pareto_k > 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

loo_compare_table <- function(loo_list, response) {
  comparison <- loo::loo_compare(x = loo_list)
  comparison <- data.frame(
    model = rownames(comparison),
    response = response,
    as.data.frame(comparison),
    row.names = NULL,
    check.names = FALSE
  )
  comparison
}

elpd_diff <- function(loo_a, loo_b) {
  if (nrow(loo_a$pointwise) != nrow(loo_b$pointwise)) {
    stop("LOO objects must have the same number of pointwise terms.")
  }

  diff_i <- loo_a$pointwise[, "elpd_loo"] - loo_b$pointwise[, "elpd_loo"]
  c(
    elpd_diff = sum(diff_i),
    se_diff = sqrt(length(diff_i) * stats::var(diff_i))
  )
}

contrast_rows <- function(loo_by_response, response) {
  pairs <- data.frame(
    hierarchical = c("ftp_sp", "ftp_sp", "pft_sp", "pft_sp"),
    one_level = c("ftp", "sp_c", "pft", "sp_c"),
    stringsAsFactors = FALSE
  )

  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    cmp <- elpd_diff(
      loo_by_response[[response]][[pairs$hierarchical[i]]],
      loo_by_response[[response]][[pairs$one_level[i]]]
    )
    data.frame(
      response = response,
      hierarchical = pairs$hierarchical[i],
      one_level = pairs$one_level[i],
      elpd_diff_hier_minus_one = cmp[["elpd_diff"]],
      se_diff = cmp[["se_diff"]],
      better = ifelse(cmp[["elpd_diff"]] > 0, "hierarchical", "one_level"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

message("Loading model fits")
fits <- lapply(model_files, readRDS)

loo_by_response <- list(
  y1 = mapply(
    FUN = function(fit, model) compute_loo(fit, model, "y1"),
    fit = fits,
    model = names(fits),
    SIMPLIFY = FALSE
  ),
  y2_observed = mapply(
    FUN = function(fit, model) {
      y2_data <- fit$data[!is.na(fit$data$y2), , drop = FALSE]
      compute_loo(fit, model, "y2", newdata = y2_data)
    },
    fit = fits,
    model = names(fits),
    SIMPLIFY = FALSE
  )
)

metrics <- do.call(
  rbind,
  unlist(
    Map(
      function(response, loos) {
        Map(
          function(model, loo_object) loo_metrics(loo_object, model, response),
          names(loos),
          loos
        )
      },
      names(loo_by_response),
      loo_by_response
    ),
    recursive = FALSE
  )
)
metrics <- merge(model_info, metrics, by = "model", all.y = TRUE)
metrics <- metrics[order(metrics$response, -metrics$elpd_loo), ]
write.csv(metrics, file.path(out_dir, "loo_metrics_by_response.csv"), row.names = FALSE)

comparisons <- do.call(
  rbind,
  Map(loo_compare_table, loo_by_response, names(loo_by_response))
)
write.csv(comparisons, file.path(out_dir, "loo_compare_by_response.csv"), row.names = FALSE)

total_metrics <- aggregate(
  cbind(elpd_loo, p_loo, looic, n_units, n_pareto_k_gt_0.7, n_pareto_k_gt_1) ~ model,
  data = metrics,
  FUN = sum
)
total_metrics$max_pareto_k <- tapply(metrics$max_pareto_k, metrics$model, max)[total_metrics$model]
total_metrics <- merge(model_info, total_metrics, by = "model", all.y = TRUE)
total_metrics$elpd_diff_from_best <- total_metrics$elpd_loo - max(total_metrics$elpd_loo)
total_metrics$looic_diff_from_best <- total_metrics$looic - min(total_metrics$looic)
total_metrics <- total_metrics[order(-total_metrics$elpd_loo), ]
write.csv(total_metrics, file.path(out_dir, "loo_total_metrics.csv"), row.names = FALSE)

contrasts <- do.call(
  rbind,
  Map(contrast_rows, list(loo_by_response), names(loo_by_response))
)
total_contrasts <- aggregate(
  cbind(elpd_diff_hier_minus_one) ~ hierarchical + one_level,
  data = contrasts,
  FUN = sum
)
total_contrasts$se_diff <- NA_real_
for (i in seq_len(nrow(total_contrasts))) {
  parts <- contrasts[
    contrasts$hierarchical == total_contrasts$hierarchical[i] &
      contrasts$one_level == total_contrasts$one_level[i],
  ]
  total_contrasts$se_diff[i] <- sqrt(sum(parts$se_diff^2))
}
total_contrasts$response <- "total_y1_plus_observed_y2"
total_contrasts$better <- ifelse(
  total_contrasts$elpd_diff_hier_minus_one > 0,
  "hierarchical",
  "one_level"
)
contrasts <- rbind(
  contrasts[, names(total_contrasts)],
  total_contrasts
)
contrasts <- contrasts[
  order(contrasts$response, contrasts$hierarchical, contrasts$one_level),
]
write.csv(contrasts, file.path(out_dir, "hierarchy_contrasts.csv"), row.names = FALSE)

saveRDS(
  list(
    model_info = model_info,
    loo_by_response = loo_by_response,
    metrics = metrics,
    comparisons = comparisons,
    total_metrics = total_metrics,
    contrasts = contrasts
  ),
  file.path(out_dir, "loo_model_comparison.rds")
)

message("\nLOO ranking by response:")
print(comparisons)

message("\nTotal ranking, summing y1 and observed y2 ELPD:")
print(total_metrics)

message("\nHierarchical versus one-level contrasts:")
print(contrasts)

message("\nWrote outputs to: ", normalizePath(out_dir))
