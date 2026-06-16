## install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop(
    "Usage: Rscript scripts/run_bef_zsh.R ",
    "<data_file> <hierarchy> <chains> <cores> <iter> <adapt_delta> <max_treedepth> <x_var> ",
    "[model] [k_hierarchy_depth] [family] [scale_x]\n",
    "Example: Rscript scripts/run_bef_zsh.R ",
    "processed_data/plot_biomass.txt 'PFT,sp_code' 4 4 4000 0.99 15 rsd\n",
    "New model example: Rscript scripts/run_bef_zsh.R ",
    "processed_data/plot_biomass.txt 'PFT,sp_code' 4 4 4000 0.99 15 rsd ",
    "log_partial_pool 1 student FALSE"
  )
}

data_file <- args[1]
hierarchy_arg <- args[2]
chains_n <- as.numeric(args[3])
cores_n <- as.numeric(args[4])
iter_n <- as.numeric(args[5])
adapt_delta_n <- as.numeric(args[6])
max_treedepth_n <- as.numeric(args[7])
x_var <- args[8]
model_type <- if (length(args) >= 9) args[9] else "exp_decay"
k_hierarchy_depth_n <- if (length(args) >= 10) as.numeric(args[10]) else 1
family_arg <- if (length(args) >= 11) args[11] else "student"
scale_x <- if (length(args) >= 12) tolower(args[12]) %in% c("1", "true", "t", "yes", "y") else FALSE

model_type <- tolower(model_type)
family_arg <- tolower(family_arg)

if (!file.exists(data_file)) {
  stop("Data file does not exist: ", data_file)
}

if (any(is.na(c(chains_n, cores_n, iter_n, adapt_delta_n, max_treedepth_n)))) {
  stop("One or more numeric arguments are NA. Check chains, cores, iter, adapt_delta, and max_treedepth.")
}

if (is.na(k_hierarchy_depth_n)) {
  stop("k_hierarchy_depth is NA.")
}

bp <- read.table(
  data_file,
  sep = "\t",
  header = TRUE
)

hierarchy_vec <- strsplit(hierarchy_arg, ",")[[1]]
hierarchy_vec <- trimws(hierarchy_vec)


library(brms)
library(cmdstanr)

# Source order matters: 00_bayesian_functions.R contains the new function,
# and 00_bayesian_exp_decay.R keeps the existing hierarchical exp_decay default.
source("scripts/00_bayesian_functions.R")
source("scripts/00_bayesian_exp_decay.R")

family_obj <- switch(
  family_arg,
  student = student(),
  gaussian = gaussian(),
  stop("Unknown family: ", family_arg, ". Use 'student' or 'gaussian'.")
)

fit_model <- function() {
  if (model_type %in% c("exp_decay", "standard", "original")) {
    return(
      Bef_bayes_exp_decay(
        bp,
        y_1 = "befa.st",
        y_2 = "befr.st",
        x = x_var,
        hierarchy = hierarchy_vec,
        chains = chains_n,
        cores = cores_n,
        iter = iter_n,
        adapt_delta = adapt_delta_n,
        max_treedepth = max_treedepth_n
      )
    )
  }

  if (model_type %in% c("log_partial_pool", "log_partial", "partial_pool")) {
    return(
      Bef_bayes_exp_decay_log_partial_pool(
        bp,
        y_1 = "befa.st",
        y_2 = "befr.st",
        x = x_var,
        hierarchy = hierarchy_vec,
        k_hierarchy_depth = k_hierarchy_depth_n,
        chains = chains_n,
        cores = cores_n,
        iter = iter_n,
        adapt_delta = adapt_delta_n,
        max_treedepth = max_treedepth_n,
        family = family_obj,
        scale_x = scale_x
      )
    )
  }

  stop(
    "Unknown model: ", model_type,
    ". Use 'exp_decay' or 'log_partial_pool'."
  )
}

tryCatch({

  start_time <- Sys.time()

  cat("\n=====================================\n")
  cat("START:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Data file:", data_file, "\n")
  cat("Rows:", nrow(bp), "\n")
  cat("Columns:", ncol(bp), "\n")
  cat("Hierarchy:", paste(hierarchy_vec, collapse = " / "), "\n")
  cat("X variable:", x_var, "\n")
  cat("Chains:", chains_n, "\n")
  cat("Iterations:", iter_n, "\n")
  cat("Cores:", cores_n, "\n")
  cat("Adapt delta:", adapt_delta_n, "\n")
  cat("Max treedepth:", max_treedepth_n, "\n")
  cat("Model:", model_type, "\n")
  cat("k hierarchy depth:", k_hierarchy_depth_n, "\n")
  cat("Family:", family_arg, "\n")
  cat("Scale x:", scale_x, "\n")
  cat("=====================================\n\n")

	fit_exp <- fit_model()
	hierarchy_name <- paste(hierarchy_vec, collapse = "_")
	hierarchy_name <- gsub("[^A-Za-z0-9_]+", "_", hierarchy_name)

	end_time <- Sys.time()
	cat("\nSUCCESS\n")
	cat("Finished:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
	cat("Elapsed:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
	cat("Hierarchy:", paste(hierarchy_vec, collapse = " / "), "\n")
	cat("X variable:", x_var, "\n")
	cat("Model:", model_type, "\n")

}, error = function(e) {

	cat("\nERROR OCCURRED\n")
	cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
	cat("Message:\n")
	cat(conditionMessage(e), "\n")
	quit(status = 1)

})
