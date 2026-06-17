bayes_none_values <- c("", "none", "null", "na", "false", "0")

`%||%` <- function(x, y) {
	if (is.null(x)) y else x
}

bayes_is_none <- function(x, include_all = FALSE) {
	if (length(x) == 0 || is.na(x)) {
		return(TRUE)
	}
	value <- tolower(trimws(as.character(x)[1]))
	none_values <- bayes_none_values
	if (include_all) {
		none_values <- c(none_values, "all")
	}
	value %in% none_values
}

bayes_parse_csv <- function(x, include_all = FALSE) {
	if (bayes_is_none(x, include_all = include_all)) {
		return(character(0))
	}
	values <- trimws(strsplit(as.character(x)[1], ",", fixed = TRUE)[[1]])
	values[nzchar(values)]
}

bayes_label <- function(values, sep = " / ", none = "none") {
	if (length(values) == 0) {
		return(none)
	}
	paste(values, collapse = sep)
}

bayes_safe_name <- function(values, none = "none") {
	if (length(values) == 0 || all(!nzchar(as.character(values)))) {
		return(none)
	}
	value <- paste(values, collapse = "-")
	value <- gsub("[^A-Za-z0-9_-]+", "-", value)
	value <- gsub("-+", "-", value)
	value <- gsub("(^-+|-+$)", "", value)
	if (!nzchar(value)) {
		value <- none
	}
	value
}

bayes_model_abbr <- function(shape) {
	if (!is.null(shape$abbr) && nzchar(shape$abbr)) {
		return(bayes_safe_name(shape$abbr))
	}
	bayes_safe_name(shape$name)
}

bayes_compact_hierarchy_name <- function(hierarchy) {
	if (length(hierarchy) == 0) {
		return("none")
	}
	paste(vapply(hierarchy, bayes_safe_name, character(1), none = "none"), collapse = "-")
}

bayes_compact_family_name <- function(family) {
	if (is.character(family)) {
		family <- family[1]
	} else if (!is.null(family$family)) {
		family <- family$family
	} else {
		family <- class(family)[1]
	}
	family <- tolower(trimws(as.character(family)[1]))
	switch(
		family,
		student = "tdis",
		stud = "tdis",
		tdis = "tdis",
		gaussian = "ndis",
		gausian = "ndis",
		gaus = "ndis",
		normal = "ndis",
		ndis = "ndis",
		bayes_safe_name(family)
	)
}

bayes_compact_depth_name <- function(k_hierarchy_depth) {
	k_hierarchy_depth <- as.integer(k_hierarchy_depth)
	if (is.na(k_hierarchy_depth)) {
		return("k0")
	}
	if (k_hierarchy_depth == 0) {
		return("k0")
	}
	paste0("k", k_hierarchy_depth)
}

bayes_compact_iter_name <- function(iter) {
	iter <- as.integer(iter)
	if (is.na(iter)) {
		return("iter")
	}
	if (iter >= 1000 && iter %% 1000 == 0) {
		return(paste0(iter / 1000, "k"))
	}
	if (iter >= 100 && iter %% 100 == 0) {
		return(paste0(iter / 100, "c"))
	}
	as.character(iter)
}

bayes_compact_chain_name <- function(chains) {
	as.character(as.integer(chains))
}

bayes_compact_delta_name <- function(adapt_delta) {
	value <- formatC(as.numeric(adapt_delta), format = "f", digits = 2)
	sub("^0\\.", "", value)
}

bayes_compact_treedepth_name <- function(max_treedepth) {
	as.character(as.integer(max_treedepth))
}

bayes_file_suffix <- function(file_suffix) {
	if (is.null(file_suffix) || bayes_is_none(file_suffix)) {
		return("")
	}
	suffix <- bayes_safe_name(file_suffix, none = "")
	if (!nzchar(suffix)) {
		return("")
	}
	paste0("_", suffix)
}

bayes_build_output_file <- function(shape,
                                    hierarchy,
                                    k_hierarchy_depth,
                                    x,
                                    family,
                                    chains,
                                    iter,
                                    adapt_delta,
                                    max_treedepth,
                                    scale_x,
                                    file_suffix = NULL) {
	paste0(
		bayes_model_abbr(shape), "_",
		bayes_compact_hierarchy_name(hierarchy), "_",
		bayes_compact_depth_name(k_hierarchy_depth), "_",
		bayes_compact_family_name(family), "_",
		bayes_safe_name(if (scale_x) paste0(x, "_scaled") else x), "_",
		bayes_compact_chain_name(chains), "-",
		bayes_compact_iter_name(iter), "-",
		bayes_compact_delta_name(adapt_delta), "-",
		bayes_compact_treedepth_name(max_treedepth),
		bayes_file_suffix(file_suffix),
		".rds"
	)
}

bayes_fit_file <- function(fit) {
	if (is.list(fit) && !is.null(fit$file)) {
		return(fit$file)
	}
	NA_character_
}

bayes_check_columns <- function(data, vars) {
	missing_vars <- setdiff(vars, names(data))
	if (length(missing_vars) > 0) {
		stop("Missing columns in data: ", paste(missing_vars, collapse = ", "))
	}
}

bayes_add_hierarchy_factors <- function(df, hierarchy) {
	for (v in hierarchy) {
		df[[v]] <- factor(df[[v]])
	}
	for (i in seq_along(hierarchy)) {
		df[[paste0("h", i)]] <- factor(
			interaction(df[, hierarchy[seq_len(i)], drop = FALSE], drop = TRUE))
	}
	df
}

bayes_random_intercept_formula <- function(depth) {
	if (depth == 0) {
		return(~ 1)
	}
	random_terms <- paste0("(1 | h", seq_len(depth), ")", collapse = " + ")
	as.formula(paste("~ 1 +", random_terms))
}

bayes_parameter_depth <- function(parameter, hierarchy_depth, k_hierarchy_depth) {
	pooling <- parameter$pooling %||% "full"
	switch(
		pooling,
		full = hierarchy_depth,
		k_depth = k_hierarchy_depth,
		none = 0,
		stop("Unknown parameter pooling type: ", pooling)
	)
}

bayes_make_prior <- function(spec, ...) {
	prior_text <- switch(
		spec$dist,
		normal = sprintf("normal(%s, %s)", signif(spec$mean, 8), signif(spec$sd, 8)),
		exponential = sprintf("exponential(%s)", signif(spec$rate, 8)),
		stop("Unknown prior distribution: ", spec$dist)
	)
	prior_string(prior_text, ...)
}

bayes_combine_priors <- function(priors) {
	if (length(priors) == 0) {
		stop("No priors were generated.")
	}
	do.call(c, priors)
}

bayes_shape_parameter_prior <- function(shape, response, parameter) {
	if (!is.null(shape$priors[[response]][[parameter]])) {
		return(shape$priors[[response]][[parameter]])
	}
	if (!is.null(shape$priors$default[[parameter]])) {
		return(shape$priors$default[[parameter]])
	}
	stop("No prior defined for response ", response, " parameter ", parameter, ".")
}

bayes_build_priors <- function(shape, responses, hierarchy_depth, k_hierarchy_depth) {
	priors <- list()
	parameter_names <- names(shape$parameters)

	for (response in responses) {
		for (parameter in parameter_names) {
			priors[[length(priors) + 1L]] <- bayes_make_prior(
				bayes_shape_parameter_prior(shape, response, parameter),
				nlpar = parameter,
				resp = response
			)
		}

		if (!is.null(shape$response_priors$sigma)) {
			priors[[length(priors) + 1L]] <- bayes_make_prior(
				shape$response_priors$sigma,
				class = "sigma",
				resp = response
			)
		}
	}

	for (parameter in parameter_names) {
		depth <- bayes_parameter_depth(shape$parameters[[parameter]], hierarchy_depth, k_hierarchy_depth)
		if (depth > 0 && !is.null(shape$parameters[[parameter]]$sd_prior)) {
			for (response in responses) {
				priors[[length(priors) + 1L]] <- bayes_make_prior(
					shape$parameters[[parameter]]$sd_prior,
					class = "sd",
					nlpar = parameter,
					resp = response
				)
			}
		}
	}

	bayes_combine_priors(priors)
}

bayes_build_nl_formula <- function(response, shape, parameter_formulas, missing = FALSE) {
	lhs <- if (missing) paste0(response, " | mi()") else response
	args <- c(
		list(formula = as.formula(paste(lhs, "~", shape$mean))),
		parameter_formulas,
		list(nl = TRUE)
	)
	do.call(bf, args)
}

Bef_bayes_fit_model_shape <- function(data,
                                      y_1,
                                      y_2,
                                      x,
                                      model_shape = NULL,
                                      hierarchy = c("ftp", "sp_code"),
                                      k_hierarchy_depth = 1,
                                      chains = 4,
                                      iter = 4000,
                                      cores = 4,
                                      adapt_delta = 0.99,
                                      max_treedepth = 15,
                                      family = Gamma(link = "identity"),
                                      scale_x = FALSE,
                                      file_refit = "on_change",
                                      file_suffix = NULL,
                                      ...) {
	if (!exists("bayes_get_model_shape", mode = "function")) {
		stop("Source scripts/00_model_shapes.R before calling Bef_bayes_fit_model_shape().")
	}
	if (is.null(model_shape)) {
		model_shape <- bayes_active_model_shape
	}
	shape <- bayes_get_model_shape(model_shape)

	hierarchy <- trimws(hierarchy)
	hierarchy <- hierarchy[nzchar(hierarchy)]

	vars <- c(y_1, y_2, x, hierarchy)
	bayes_check_columns(data, vars)

	k_hierarchy_depth <- as.integer(k_hierarchy_depth)
	if (length(k_hierarchy_depth) != 1 || is.na(k_hierarchy_depth)) {
		stop("k_hierarchy_depth must be one integer.")
	}
	if (k_hierarchy_depth < 0 || k_hierarchy_depth > length(hierarchy)) {
		stop("k_hierarchy_depth must be between 0 and length(hierarchy).")
	}

	df <- data[, vars, drop = FALSE]
	names(df) <- c("y1", "y2", "x", hierarchy)
	df <- df[complete.cases(df[, c("y1", "x", hierarchy), drop = FALSE]), , drop = FALSE]
	if (nrow(df) == 0) {
		stop("No complete rows for y1, x, and hierarchy columns.")
	}

	if (any(df$y1 <= 1, na.rm = TRUE)) {
		stop(y_1, " must be > 1 for log(y1 - 1).")
	}
	if (any(!is.na(df$y2) & df$y2 <= 0)) {
		stop(y_2, " must be > 0 for log(y2).")
	}

	df$y1m1 <- df$y1 - 1
	if (any(df$y1m1 <= 0, na.rm = TRUE)) {
		stop(y_1, " - 1 must be > 0 for original-scale modeling.")
	}
	df <- bayes_add_hierarchy_factors(df, hierarchy)

	if (scale_x) {
		x_scale <- median(df$x, na.rm = TRUE)
		if (!is.finite(x_scale) || x_scale == 0) {
			stop("Cannot scale x because its median is not finite and non-zero.")
		}
		df$x <- df$x / x_scale
	} else {
		x_scale <- 1
	}

	parameter_formulas <- lapply(shape$parameters, function(parameter) {
		depth <- bayes_parameter_depth(parameter, length(hierarchy), k_hierarchy_depth)
		bayes_random_intercept_formula(depth)
	})

	bf_y1 <- bayes_build_nl_formula("y1m1", shape, parameter_formulas)
	bf_y2 <- bayes_build_nl_formula("y2", shape, parameter_formulas, missing = TRUE)
	priors <- bayes_build_priors(
		shape,
		responses = c("y1m1", "y2"),
		hierarchy_depth = length(hierarchy),
		k_hierarchy_depth = k_hierarchy_depth
	)

	dir.create("bayes_outputs", recursive = TRUE, showWarnings = FALSE)
	output_file <- file.path(
		"bayes_outputs",
		bayes_build_output_file(
			shape = shape,
			hierarchy = hierarchy,
			k_hierarchy_depth = k_hierarchy_depth,
			x = x,
			family = family,
			chains = chains,
			iter = iter,
			adapt_delta = adapt_delta,
			max_treedepth = max_treedepth,
			scale_x = scale_x,
			file_suffix = file_suffix
		)
	)

	fit <- brm(
		bf_y1 +
			bf_y2 + set_rescor(FALSE),
		data = df,
		family = family,
		prior = priors,
		chains = chains,
		iter = iter,
		cores = cores,
		control = list(
			adapt_delta = adapt_delta,
			max_treedepth = max_treedepth),
		backend = "cmdstanr",
		file = output_file,
		file_refit = file_refit,
		save_pars = save_pars(all = TRUE),
		...
	)

	attr(fit, "x_scale") <- x_scale
	attr(fit, "k_hierarchy_depth") <- k_hierarchy_depth
	attr(fit, "model_shape") <- shape$name
	return(fit)
}

Bef_bayes_exp_decay_log_partial_pool <- function(...) {
	Bef_bayes_fit_model_shape(..., model_shape = "exp_decay")
}
