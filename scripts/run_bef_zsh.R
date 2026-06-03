## install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop(
    "Usage: Rscript scripts/run_bef_zsh.R ",
    "<data_file> <hierarchy> <chains> <cores> <iter> <adapt_delta> <max_treedepth> <x_var>\n",
    "Example: Rscript scripts/run_bef_zsh.R ",
    "processed_data/plot_biomass.txt 'PFT,sp_code' 4 4 4000 0.99 15 rsd"
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

if (!file.exists(data_file)) {
  stop("Data file does not exist: ", data_file)
}

if (any(is.na(c(chains_n, cores_n, iter_n, adapt_delta_n, max_treedepth_n)))) {
  stop("One or more numeric arguments are NA. Check chains, cores, iter, adapt_delta, and max_treedepth.")
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

source("scripts/00_bayesian_exp_decay.R")

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
  cat("=====================================\n\n")

	fit_exp <- Bef_bayes_exp_decay(
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
	hierarchy_name <- paste(hierarchy_vec, collapse = "_")
	hierarchy_name <- gsub("[^A-Za-z0-9_]+", "_", hierarchy_name)

	end_time <- Sys.time()
	cat("\nSUCCESS\n")
	cat("Finished:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
	cat("Elapsed:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
	cat("Hierarchy:", paste(hierarchy_vec, collapse = " / "), "\n")
	cat("X variable:", x_var, "\n")

}, error = function(e) {

	cat("\nERROR OCCURRED\n")
	cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
	cat("Message:\n")
	cat(conditionMessage(e), "\n")
	quit(status = 1)

})
