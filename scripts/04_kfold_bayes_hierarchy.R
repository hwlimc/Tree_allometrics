## Exact K-fold cross-validation for the Bayesian hierarchy models.
##
## Recommended smoke test:
## K=2 KSUB=1 MODELS=ftp_sp CHAINS=1 ITER=1000 CORES=1 \
##   Rscript --vanilla scripts/04_kfold_bayes_hierarchy.R
##
## Full comparison:
## K=10 MODELS=ftp_sp,pft_sp,sp_c,pft,ftp CHAINS=4 ITER=4000 CORES=4 \
##   Rscript --vanilla scripts/04_kfold_bayes_hierarchy.R

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

parse_bool <- function(x, default = FALSE) {
  if (!nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_int_vec <- function(x, default = NULL) {
  if (!nzchar(x)) {
    return(default)
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.integer(value) else default
}

K <- env_int("K", 10L)
Ksub <- parse_int_vec(Sys.getenv("KSUB", unset = ""), default = seq_len(K))
seed <- env_int("SEED", 20260616L)
chains <- env_int("CHAINS", 4L)
iter <- env_int("ITER", 4000L)
cores <- env_int("CORES", 4L)
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", unset = "0.99"))
max_treedepth <- env_int("MAX_TREEDEPTH", 15L)
refresh <- env_int("REFRESH", 250L)
force_refit <- parse_bool(Sys.getenv("FORCE_REFIT", unset = ""), default = FALSE)
save_fold_fits <- parse_bool(Sys.getenv("SAVE_FOLD_FITS", unset = ""), default = FALSE)
recompile <- parse_bool(Sys.getenv("RECOMPILE", unset = ""), default = FALSE)
sample_new_levels <- Sys.getenv("SAMPLE_NEW_LEVELS", unset = "gaussian")
fold_mode <- Sys.getenv("FOLD_MODE", unset = "stratified_y1_y2")

models_requested <- Sys.getenv(
  "MODELS",
  unset = "ftp_sp,pft_sp,sp_c,pft,ftp"
)
models_to_run <- strsplit(models_requested, ",", fixed = TRUE)[[1]]
models_to_run <- trimws(models_to_run)
unknown_models <- setdiff(models_to_run, names(model_files))
if (length(unknown_models) > 0) {
  stop("Unknown MODELS entries: ", paste(unknown_models, collapse = ", "))
}

out_dir <- Sys.getenv("OUT_DIR", unset = "bayes_outputs/kfold_cv")
score_dir <- file.path(out_dir, "fold_scores")
fit_dir <- file.path(out_dir, "fold_fits")
dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
if (save_fold_fits) {
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
}

log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(mean(exp(x - max_x)))
}

score_response <- function(fit, newdata, response, sample_new_levels) {
  if (nrow(newdata) == 0) {
    return(numeric(0))
  }

  log_lik_matrix <- brms::log_lik(
    fit,
    newdata = newdata,
    resp = response,
    sample_new_levels = sample_new_levels
  )
  if (anyNA(log_lik_matrix)) {
    stop("NA values found in held-out log likelihood for ", response)
  }
  apply(log_lik_matrix, 2, log_mean_exp)
}

kfold_summary <- function(pointwise, response) {
  x <- pointwise$elpd[pointwise$response == response]
  data.frame(
    response = response,
    folds_evaluated = length(unique(pointwise$fold)),
    n_units = length(x),
    elpd_kfold = sum(x),
    se_elpd_kfold = sqrt(length(x) * stats::var(x)),
    stringsAsFactors = FALSE
  )
}

make_folds <- function(data, K, seed, mode) {
  set.seed(seed)
  if (identical(mode, "random")) {
    return(loo::kfold_split_random(K = K, N = nrow(data)))
  }
  if (identical(mode, "stratified_y2")) {
    return(loo::kfold_split_stratified(K = K, x = !is.na(data$y2)))
  }
  if (identical(mode, "stratified_y1_y2")) {
    probs <- seq(0, 1, length.out = min(6L, max(3L, K)) + 1L)
    breaks <- unique(as.numeric(stats::quantile(data$y1, probs, na.rm = TRUE)))
    y1_bin <- cut(data$y1, breaks = breaks, include.lowest = TRUE)
    strata <- interaction(y1_bin, !is.na(data$y2), drop = TRUE)
    return(loo::kfold_split_stratified(K = K, x = strata))
  }
  stop("Unknown FOLD_MODE: ", mode)
}

score_file <- function(model, fold) {
  file.path(score_dir, paste0("kfold_", model, "_fold_", fold, ".rds"))
}

fold_fit_file <- function(model, fold) {
  if (!save_fold_fits) {
    return(NULL)
  }
  file.path(fit_dir, paste0("kfold_", model, "_fold_", fold))
}

fit_one_fold <- function(base_fit, model, fold, folds) {
  cache_file <- score_file(model, fold)
  if (!force_refit && file.exists(cache_file)) {
    message("Using cached score: ", model, " / fold ", fold)
    return(readRDS(cache_file))
  }

  message("Fitting ", model, " / fold ", fold, " of ", K)
  train_data <- base_fit$data[folds != fold, , drop = FALSE]
  test_data <- base_fit$data[folds == fold, , drop = FALSE]

  set.seed(seed + fold)
  fold_fit <- update(
    base_fit,
    newdata = train_data,
    recompile = recompile,
    chains = chains,
    iter = iter,
    cores = cores,
    seed = seed + fold,
    refresh = refresh,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    file = fold_fit_file(model, fold),
    file_refit = "on_change",
    save_pars = save_pars(all = TRUE)
  )

  y1_elpd <- score_response(
    fold_fit,
    newdata = test_data,
    response = "y1",
    sample_new_levels = sample_new_levels
  )

  y2_test_data <- test_data[!is.na(test_data$y2), , drop = FALSE]
  y2_elpd <- score_response(
    fold_fit,
    newdata = y2_test_data,
    response = "y2",
    sample_new_levels = sample_new_levels
  )

  pointwise <- rbind(
    data.frame(
      model = model,
      fold = fold,
      response = "y1",
      row_id = which(folds == fold),
      elpd = y1_elpd,
      stringsAsFactors = FALSE
    ),
    data.frame(
      model = model,
      fold = fold,
      response = "y2_observed",
      row_id = which(folds == fold)[!is.na(test_data$y2)],
      elpd = y2_elpd,
      stringsAsFactors = FALSE
    )
  )

  out <- list(
    model = model,
    fold = fold,
    train_n = nrow(train_data),
    test_n = nrow(test_data),
    test_y2_observed_n = nrow(y2_test_data),
    pointwise = pointwise
  )
  saveRDS(out, cache_file)
  out
}

message("Loading reference fit for fold construction")
reference_fit <- readRDS(model_files[[models_to_run[1]]])
folds <- make_folds(reference_fit$data, K = K, seed = seed, mode = fold_mode)
fold_table <- data.frame(row_id = seq_along(folds), fold = folds)
write.csv(fold_table, file.path(out_dir, "fold_assignments.csv"), row.names = FALSE)

message("K = ", K)
message("Folds to evaluate = ", paste(Ksub, collapse = ", "))
message("Fold mode = ", fold_mode)
message("Models = ", paste(models_to_run, collapse = ", "))
message("chains = ", chains, ", iter = ", iter, ", cores = ", cores)

all_scores <- list()
for (model in models_to_run) {
  message("\nLoading base fit: ", model)
  base_fit <- readRDS(model_files[[model]])
  if (nrow(base_fit$data) != length(folds)) {
    stop("Model ", model, " has a different number of rows than the fold vector.")
  }
  fold_scores <- lapply(Ksub, fit_one_fold, base_fit = base_fit, model = model, folds = folds)
  all_scores[[model]] <- do.call(rbind, lapply(fold_scores, `[[`, "pointwise"))
}

pointwise <- do.call(rbind, all_scores)
write.csv(pointwise, file.path(out_dir, "kfold_pointwise.csv"), row.names = FALSE)

by_response <- do.call(
  rbind,
  lapply(split(pointwise, pointwise$model), function(df) {
    do.call(rbind, lapply(c("y1", "y2_observed"), kfold_summary, pointwise = df))
  })
)
by_response$model <- rep(names(split(pointwise, pointwise$model)), each = 2L)
by_response <- merge(model_info, by_response, by = "model", all.y = TRUE)
by_response <- by_response[order(by_response$response, -by_response$elpd_kfold), ]
write.csv(by_response, file.path(out_dir, "kfold_metrics_by_response.csv"), row.names = FALSE)

total <- aggregate(
  cbind(elpd_kfold, n_units) ~ model,
  data = by_response,
  FUN = sum
)
total$se_elpd_kfold <- NA_real_
for (i in seq_len(nrow(total))) {
  model_rows <- by_response[by_response$model == total$model[i], ]
  total$se_elpd_kfold[i] <- sqrt(sum(model_rows$se_elpd_kfold^2))
}
total <- merge(model_info, total, by = "model", all.y = TRUE)
total$elpd_diff_from_best <- total$elpd_kfold - max(total$elpd_kfold)
total <- total[order(-total$elpd_kfold), ]
write.csv(total, file.path(out_dir, "kfold_total_metrics.csv"), row.names = FALSE)

message("\nK-fold metrics by response:")
print(by_response)

message("\nK-fold total metrics:")
print(total)

message("\nWrote outputs to: ", normalizePath(out_dir))
