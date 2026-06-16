### Bayesian parameterization following a Michaelis-Menten function
# # Bef_bayes_mich_men<- function(data,
                                # # y_1,
                                # # y_2,
                                # # x,
                                # # group,
                                # # chains = 4,
                                # # iter = 4000,
                                # # cores = 4,
                                # # adapt_delta = 0.95,
                                # # max_treedepth = 12) {
	# # df <- data[, c(y_1, y_2, x, group)]
	# # names(df) <- c("y1", "y2", "x", "grp")

	# # df <- df[complete.cases(df[, c("x", "grp")]), ]
	# # df$grp <- factor(df$grp)
	# # med_x <- median(df$x, na.rm = TRUE)
	# # df$x <- df$x / med_x
	# # prior_r0 <- "lognormal(0, 1)"
	
	# # bf_y1 <- bf(
		# # y1 ~ L + A / (1 + x / r0),
		# # L  ~ 1 + (1 | grp),
		# # A  ~ 1 + (1 | grp),
		# # r0 ~ 1 + (1 | grp),
		# # nl = TRUE)

	# # bf_y2 <- bf(
		# # y2 | mi() ~ L + A / (1 + x / r0),
		# # L  ~ 1 + (1 | grp),
		# # A  ~ 1 + (1 | grp),
		# # r0 ~ 1 + (1 | grp),
		# # nl = TRUE)

# # fit <- brm(
	# # bf_y1 + bf_y2 + set_rescor(FALSE),
	# # data = df,
	# # family = gaussian(),
	# # prior = c(
		# # prior(normal(1.2, 0.5), nlpar = "L",  resp = "y1", lb = 0),
		# # prior(normal(0.6, 0.6), nlpar = "A",  resp = "y1", lb = 0),
		# # prior_string(prior_r0, nlpar = "r0", resp = "y1", lb = 0),

		# # prior(normal(1.4, 0.5), nlpar = "L",  resp = "y2", lb = 0),
		# # prior(normal(0.7, 0.7), nlpar = "A",  resp = "y2", lb = 0),
		# # prior_string(prior_r0, nlpar = "r0", resp = "y2", lb = 0)),
	# # chains = chains,
	# # iter = iter,
	# # cores = cores,
	# # control = list(
		# # adapt_delta = adapt_delta,
		# # max_treedepth = max_treedepth),
	# # backend = "cmdstanr",
	# # file = paste0("processed_data/Bayes_mich_men_",
	# # group,"_",
	# # x,"_",
    # # chains, "chn_",
    # # iter, "itr_",
	# # cores, "cor_", 
    # # adapt_delta, "del_",
    # # max_treedepth, "depth",
    # # ".rds"),
	# # file_refit = "on_change",
	# # save_pars = save_pars(all = TRUE)
	# # )
	# # return(fit)
# # }


### Bayesian parameterization following an exponential decay function

Bef_bayes_exp_decay<- function(data,
                                y_1,
                                y_2,
                                x,
                                group,
                                chains = 4,
                                iter = 4000,
                                cores = 4,
                                adapt_delta = 0.99,
                                max_treedepth = 15) {
	df <- data[, c(y_1, y_2, x, group)]
	names(df) <- c("y1", "y2", "x", "grp")

	df <- df[complete.cases(df[, c("y1", "x", "grp")]), ]
	df$grp <- factor(df$grp)
	med_x <- median(df$x, na.rm = TRUE)
	df$x <- df$x / med_x
	prior_k <- "lognormal(0,1)"
  
	bf_y1 <- bf(
		y1 ~ L + A * exp(-k * x),
		L  ~ 1 + (1 | grp),
		A  ~ 1 + (1 | grp),
		k ~ 1 + (1 | grp),
		nl = TRUE)

	bf_y2 <- bf(
		y2 | mi() ~ L + A * exp(-k * x),
		L  ~ 1 + (1 | grp),
		A  ~ 1 + (1 | grp),
		k ~ 1 + (1 | grp),
		nl = TRUE)

fit <- brm(
	bf_y1 + bf_y2 + set_rescor(FALSE),
	data = df,
	family = gaussian(),
	prior = c(
		prior(normal(1.2, 0.5), nlpar = "L",  resp = "y1", lb = 0),
		prior(normal(0.6, 0.6), nlpar = "A",  resp = "y1", lb = 0),
		prior_string(prior_k, nlpar = "k", resp = "y1", lb = 0),

		prior(normal(1.4, 0.5), nlpar = "L",  resp = "y2", lb = 0),
		prior(normal(0.7, 0.7), nlpar = "A",  resp = "y2", lb = 0),
		prior_string(prior_k, nlpar = "k", resp = "y2", lb = 0)),
	chains = chains,
	iter = iter,
	cores = cores,
	control = list(
		adapt_delta = adapt_delta,
		max_treedepth = max_treedepth),
	backend = "cmdstanr",
	file = paste0("processed_data/Bayes_exp_decay_",
	group,"_",
	x,"_",
    chains, "chn_",
    iter, "itr_",
	cores, "cor_", 
    adapt_delta, "del_",
    max_treedepth, "depth",
    ".rds"),
	file_refit = "on_change",
	save_pars = save_pars(all = TRUE)
	)
	return(fit)
}


### Bayesian exponential decay with log parameters and simpler k hierarchy
### This is intended for sparse species-level data where species-specific k
### is weakly identified. L and A can vary across the full hierarchy, while k
### defaults to only the first hierarchy level.

Bef_bayes_exp_decay_log_partial_pool <- function(data,
                                                 y_1,
                                                 y_2,
                                                 x,
                                                 hierarchy = c("ftp", "sp_code"),
                                                 k_hierarchy_depth = 1,
                                                 chains = 4,
                                                 iter = 4000,
                                                 cores = 4,
                                                 adapt_delta = 0.99,
                                                 max_treedepth = 15,
                                                 family = student(),
                                                 scale_x = FALSE,
                                                 file_refit = "on_change",
                                                 ...) {
	vars <- c(y_1, y_2, x, hierarchy)
	missing_vars <- setdiff(vars, names(data))
	if (length(missing_vars) > 0) {
		stop("Missing columns in data: ", paste(missing_vars, collapse = ", "))
	}

	k_hierarchy_depth <- as.integer(k_hierarchy_depth)
	if (length(k_hierarchy_depth) != 1 || is.na(k_hierarchy_depth)) {
		stop("k_hierarchy_depth must be one integer.")
	}
	if (k_hierarchy_depth < 0 || k_hierarchy_depth > length(hierarchy)) {
		stop("k_hierarchy_depth must be between 0 and length(hierarchy).")
	}

	df <- data[, vars]
	names(df) <- c("y1", "y2", "x", hierarchy)
	df <- df[complete.cases(df[, c("y1", "x", hierarchy)]), ]

	for (v in hierarchy) {
		df[[v]] <- factor(df[[v]])
	}
	for (i in seq_along(hierarchy)) {
		new_name <- paste0("h", i)
		df[[new_name]] <- factor(
			interaction(df[, hierarchy[seq_len(i)], drop = FALSE], drop = TRUE))
	}

	if (scale_x) {
		x_scale <- median(df$x, na.rm = TRUE)
		df$x <- df$x / x_scale
	} else {
		x_scale <- 1
	}

	random_formula <- function(depth) {
		if (depth == 0) {
			return(~ 1)
		}
		random_terms <- paste0("(1 | h", seq_len(depth), ")", collapse = " + ")
		as.formula(paste("~ 1 +", random_terms))
	}

	LA_formula <- random_formula(length(hierarchy))
	k_formula <- random_formula(k_hierarchy_depth)

	bf_y1 <- bf(
		y1 ~ exp(logL) + exp(logA) * exp(-exp(logk) * x),
		logL = LA_formula,
		logA = LA_formula,
		logk = k_formula,
		nl = TRUE)

	bf_y2 <- bf(
		y2 | mi() ~ exp(logL) + exp(logA) * exp(-exp(logk) * x),
		logL = LA_formula,
		logA = LA_formula,
		logk = k_formula,
		nl = TRUE)

	priors <- c(
		prior(normal(log(1.2), 0.4), nlpar = "logL", resp = "y1"),
		prior(normal(log(0.6), 0.7), nlpar = "logA", resp = "y1"),
		prior(normal(log(1.0), 0.8), nlpar = "logk", resp = "y1"),

		prior(normal(log(0.75), 0.5), nlpar = "logL", resp = "y2"),
		prior(normal(log(0.7), 0.7), nlpar = "logA", resp = "y2"),
		prior(normal(log(1.0), 0.8), nlpar = "logk", resp = "y2"),

		prior(exponential(1), class = "sigma", resp = "y1"),
		prior(exponential(1), class = "sigma", resp = "y2")
	)

	if (length(hierarchy) > 0) {
		priors <- c(
			priors,
			prior(exponential(1), class = "sd", nlpar = "logL", resp = "y1"),
			prior(exponential(1), class = "sd", nlpar = "logA", resp = "y1"),
			prior(exponential(1), class = "sd", nlpar = "logL", resp = "y2"),
			prior(exponential(1), class = "sd", nlpar = "logA", resp = "y2")
		)
	}

	if (k_hierarchy_depth > 0) {
		priors <- c(
			priors,
			prior(exponential(1), class = "sd", nlpar = "logk", resp = "y1"),
			prior(exponential(1), class = "sd", nlpar = "logk", resp = "y2")
		)
	}

	hierarchy_name <- paste(hierarchy, collapse = "_")
	hierarchy_name <- gsub("[^A-Za-z0-9_]+", "_", hierarchy_name)
	family_name <- if (is.character(family)) family[1] else family$family
	family_name <- gsub("[^A-Za-z0-9_]+", "_", family_name)
	x_label <- if (scale_x) paste0(x, "_scaled") else x

	if (!dir.exists("bayes_outputs")) {
		dir.create("bayes_outputs", recursive = TRUE)
	}

	output_file <- paste0(
		"bayes_outputs/exp_decay_log_partial_pool_",
		hierarchy_name, "_",
		"kdepth", k_hierarchy_depth, "_",
		x_label, "_",
		family_name, "_",
		chains, "chn_",
		iter, "itr_",
		cores, "cor_",
		adapt_delta, "del_",
		max_treedepth, "depth",
		".rds"
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
	return(fit)
}
