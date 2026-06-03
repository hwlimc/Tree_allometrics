Bef_bayes_exp_decay <- function(data,
                                y_1,
                                y_2,
                                x,
                                hierarchy = c("ft1.forest_type", "PFT", "sp_code"),
                                chains = 4,
                                cores = 4,
								iter = 4000,
								adapt_delta = 0.95,
                                max_treedepth = 12){
	vars <- c(y_1, y_2, x, hierarchy)
	missing_vars <- setdiff(vars, names(data))
		if (length(missing_vars) > 0) {
			stop("Missing columns in data: ", paste(missing_vars, collapse = ", "))}

	df <- data[, vars]
	names(df) <- c("y1", "y2", "x", hierarchy)
	df <- df[complete.cases(df[, c("y1","x", hierarchy)]), ]

		for (v in hierarchy) {
			df[[v]] <- factor(df[[v]])
			}
		for (i in seq_along(hierarchy)) {
			new_name <- paste0("h", i)
			df[[new_name]] <- factor(
			interaction(df[, hierarchy[seq_len(i)], drop = FALSE], drop = TRUE))
			}

	nested_terms <- paste0("(1 | h", seq_along(hierarchy), ")", collapse = " + ")

	L_formula <- as.formula(paste("~ 1 +", nested_terms))
	A_formula <- as.formula(paste("~ 1 +", nested_terms))
	k_formula <- as.formula(paste("~ 1 +", nested_terms))

	prior_k <- "lognormal(0,1)"

	bf_y1 <- bf(y1 ~ L + A * exp(-k * x),
				L = L_formula,
				A = A_formula,
				k = k_formula,
				nl = TRUE
				)

	bf_y2 <- bf(y2 | mi() ~ L + A * exp(-k * x),
				L = L_formula,
				A = A_formula,
				k = k_formula,
				nl = TRUE
				)

	hierarchy_name <- paste(hierarchy, collapse = "_")
	hierarchy_name <- gsub("[^A-Za-z0-9_]+", "_", hierarchy_name)

if (!dir.exists("bayes_outputs")) {
  dir.create("bayes_outputs", recursive = TRUE)
}

	fit <- brm(
		bf_y1 + 
		bf_y2 + set_rescor(FALSE),
		data = df,
		family = gaussian(),
		prior = c(
			prior(normal(1.2, 0.5), nlpar = "L", resp = "y1", lb = 1),
			prior(normal(0.6, 0.6), nlpar = "A", resp = "y1", lb = 0),
			prior_string(prior_k, nlpar = "k", resp = "y1", lb = 0),

			prior(normal(0.75, 0.5), nlpar = "L", resp = "y2", lb = 0),
			prior(normal(0.7, 0.7), nlpar = "A", resp = "y2", lb = 0),
			prior_string(prior_k, nlpar = "k", resp = "y2", lb = 0)
				),
		chains = chains,
		iter = iter,
		cores = cores,
		control = list(
			adapt_delta = adapt_delta,
			max_treedepth = max_treedepth
			),
		backend = "cmdstanr",
		file = paste0(
			"bayes_outputs/exp_decay_",
			hierarchy_name, "_",
			x, "_",
			chains, "chn_",
			iter, "itr_",
			cores, "cor_",
			adapt_delta, "del_",
			max_treedepth, "depth",
			".rds"
			),
		file_refit = "on_change",
		save_pars = save_pars(all = TRUE)
		)

	return(fit)
	}