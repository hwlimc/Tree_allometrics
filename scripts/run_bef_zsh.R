usage <- paste0(
	"Usage: Rscript scripts/run_bef_zsh.R ",
	"<data_file> <hierarchy> <chains> <cores> <iter> <adapt_delta> <max_treedepth> <x_var> ",
	"[model] [k_hierarchy_depth] [family] [scale_x] [split_col] [split_values] [drop_split_from_hierarchy]\n",
	"Model is a shape name from scripts/00_model_shapes.R, such as exp_decay, michaelis_menten, or linear.\n",
	"Family can be gamma, lognormal/lnorm, student/tdis, or gaussian/normal/ndis.\n",
	"Combined example: Rscript scripts/run_bef_zsh.R ",
	"processed_data/plot_biomass.txt 'PFT,sp_code' 4 4 4000 0.99 15 rsd ",
	"exp_decay 1 gamma FALSE none all TRUE\n",
	"Separate species example: Rscript scripts/run_bef_zsh.R ",
	"processed_data/plot_biomass.txt none 4 4 4000 0.99 15 rsd ",
	"exp_decay 0 gamma FALSE sp_code all TRUE"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
	stop(usage)
}

arg <- function(i, default = NULL) {
	if (length(args) >= i && nzchar(args[i])) {
		return(args[i])
	}
	default
}

as_number <- function(value, name) {
	value <- as.numeric(value)
	if (is.na(value)) {
		stop(name, " is NA or not numeric.")
	}
	value
}

as_flag <- function(value, default = FALSE) {
	if (is.null(value) || is.na(value) || !nzchar(value)) {
		return(default)
	}
	tolower(trimws(value)) %in% c("1", "true", "t", "yes", "y")
}

print_fields <- function(fields) {
	for (name in names(fields)) {
		cat(name, ": ", fields[[name]], "\n", sep = "")
	}
}

data_file <- arg(1)
hierarchy_arg <- arg(2)
chains_n <- as_number(arg(3), "chains")
cores_n <- as_number(arg(4), "cores")
iter_n <- as_number(arg(5), "iter")
adapt_delta_n <- as_number(arg(6), "adapt_delta")
max_treedepth_n <- as_number(arg(7), "max_treedepth")
x_var <- arg(8)
model_type <- tolower(arg(9, "exp_decay"))
k_hierarchy_depth_n <- as_number(arg(10, "1"), "k_hierarchy_depth")
family_arg <- tolower(arg(11, "gamma"))
scale_x <- as_flag(arg(12, "FALSE"))
split_col <- trimws(arg(13, "none"))
split_values_arg <- trimws(arg(14, "all"))
drop_split_from_hierarchy <- as_flag(arg(15, "TRUE"), default = TRUE)

if (!file.exists(data_file)) {
	stop("Data file does not exist: ", data_file)
}

suppressPackageStartupMessages({
	library(brms)
	library(cmdstanr)
})
source("scripts/00_model_shapes.R")
source("scripts/00_bayesian_functions.R")
model_shape <- bayes_get_model_shape(model_type)
model_type <- model_shape$name

bp <- read.table(data_file, sep = "\t", header = TRUE)
hierarchy_vec <- bayes_parse_csv(hierarchy_arg)

has_split <- !bayes_is_none(split_col, include_all = TRUE)
if (has_split && !split_col %in% names(bp)) {
	stop("split_col is not in data: ", split_col)
}

split_values <- character(0)
if (has_split) {
	split_values <- bayes_parse_csv(split_values_arg, include_all = TRUE)
	if (length(split_values) == 0) {
		split_values <- sort(unique(as.character(bp[[split_col]])))
		split_values <- split_values[!is.na(split_values)]
	}
}

family_obj <- switch(
	family_arg,
	student = student(),
	stud = student(),
	tdis = student(),
	gamma = Gamma(link = "identity"),
	lognormal = lognormal(),
	lnorm = lognormal(),
	gaussian = gaussian(),
	gausian = gaussian(),
	gaus = gaussian(),
	normal = gaussian(),
	ndis = gaussian(),
	stop("Unknown family: ", family_arg, ". Use 'lognormal'/'lnorm', 'gamma', 'student'/'tdis', or 'gaussian'/'normal'/'ndis'.")
)

fit_model <- function(data, hierarchy, file_suffix = NULL) {
	Bef_bayes_fit_model_shape(
		data,
		y_1 = "befa.st",
		y_2 = "befr.st",
		x = x_var,
		model_shape = model_shape,
		hierarchy = hierarchy,
		k_hierarchy_depth = k_hierarchy_depth_n,
		chains = chains_n,
		cores = cores_n,
		iter = iter_n,
		adapt_delta = adapt_delta_n,
		max_treedepth = max_treedepth_n,
		family = family_obj,
		scale_x = scale_x,
		file_suffix = file_suffix
	)
}

fit_subset <- function(split_value) {
	current_data <- bp[as.character(bp[[split_col]]) == split_value, , drop = FALSE]
	if (nrow(current_data) == 0) {
		stop("No rows found for ", split_col, " = ", split_value)
	}

	current_hierarchy <- hierarchy_vec
	if (drop_split_from_hierarchy) {
		current_hierarchy <- setdiff(current_hierarchy, split_col)
	}

	cat("\n-------------------------------------\n")
	print_fields(c(
		"Subset" = paste(split_col, "=", split_value),
		"Subset rows" = nrow(current_data),
		"Subset hierarchy" = bayes_label(current_hierarchy)
	))

	fit <- fit_model(
		current_data,
		current_hierarchy,
		file_suffix = paste0(tolower(split_col), "-", split_value)
	)
	output_file <- bayes_fit_file(fit)
	rm(fit)
	invisible(gc())
	output_file
}

tryCatch({
	start_time <- Sys.time()

	cat("\n=====================================\n")
	print_fields(c(
		"START" = format(start_time, "%Y-%m-%d %H:%M:%S"),
		"Data file" = data_file,
		"Rows" = nrow(bp),
		"Columns" = ncol(bp),
		"Hierarchy" = bayes_label(hierarchy_vec),
		"X variable" = x_var,
		"Chains" = chains_n,
		"Iterations" = iter_n,
		"Cores" = cores_n,
		"Adapt delta" = adapt_delta_n,
		"Max treedepth" = max_treedepth_n,
		"Model" = model_type,
		"k hierarchy depth" = k_hierarchy_depth_n,
		"Family" = family_arg,
		"Scale x" = scale_x,
		"Split column" = if (has_split) split_col else "none"
	))
	if (has_split) {
		print_fields(c(
			"Split values" = paste(split_values, collapse = ", "),
			"Drop split from hierarchy" = drop_split_from_hierarchy
		))
	}
	cat("=====================================\n\n")

	if (has_split) {
		output_files <- setNames(rep(NA_character_, length(split_values)), split_values)
		for (split_value in split_values) {
			output_files[[split_value]] <- fit_subset(split_value)
			cat("Saved model file:", output_files[[split_value]], "\n")
		}
	} else {
		fit <- fit_model(bp, hierarchy_vec)
		output_files <- bayes_fit_file(fit)
		rm(fit)
		invisible(gc())
		cat("Saved model file:", output_files, "\n")
	}

	end_time <- Sys.time()
	cat("\nSUCCESS\n")
	print_fields(c(
		"Finished" = format(end_time, "%Y-%m-%d %H:%M:%S"),
		"Elapsed" = paste(round(difftime(end_time, start_time, units = "mins"), 2), "minutes"),
		"Hierarchy" = bayes_label(hierarchy_vec),
		"X variable" = x_var,
		"Model" = model_type
	))
	if (has_split) {
		cat("Saved model files:\n")
		for (name in names(output_files)) {
			cat("  ", name, ": ", output_files[[name]], "\n", sep = "")
		}
	} else {
		cat("Saved model file:", output_files, "\n")
	}
}, error = function(e) {
	cat("\nERROR OCCURRED\n")
	cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
	cat("Message:\n")
	cat(conditionMessage(e), "\n")
	quit(status = 1)
})
