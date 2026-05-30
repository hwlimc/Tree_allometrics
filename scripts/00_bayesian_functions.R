### Bayesian parameterization following a Michaelis-Menten function
Bef_bayes_mich_men<- function(data,
                                y_1,
                                y_2,
                                x,
                                group,
                                chains = 4,
                                iter = 4000,
                                cores = 4,
                                adapt_delta = 0.95,
                                max_treedepth = 12) {
	df <- data[, c(y_1, y_2, x, group)]
	names(df) <- c("y1", "y2", "x", "grp")

	df <- df[complete.cases(df[, c("x", "grp")]), ]
	df$grp <- factor(df$grp)
	med_x <- median(df$x, na.rm = TRUE)
	df$x <- df$x / med_x
	prior_r0 <- "lognormal(0, 1)"
	
	bf_y1 <- bf(
		y1 ~ L + A / (1 + x / r0),
		L  ~ 1 + (1 | grp),
		A  ~ 1 + (1 | grp),
		r0 ~ 1 + (1 | grp),
		nl = TRUE)

	bf_y2 <- bf(
		y2 | mi() ~ L + A / (1 + x / r0),
		L  ~ 1 + (1 | grp),
		A  ~ 1 + (1 | grp),
		r0 ~ 1 + (1 | grp),
		nl = TRUE)

fit <- brm(
	bf_y1 + bf_y2 + set_rescor(FALSE),
	data = df,
	family = gaussian(),
	prior = c(
		prior(normal(1.2, 0.5), nlpar = "L",  resp = "y1", lb = 0),
		prior(normal(0.6, 0.6), nlpar = "A",  resp = "y1", lb = 0),
		prior_string(prior_r0, nlpar = "r0", resp = "y1", lb = 0),

		prior(normal(1.4, 0.5), nlpar = "L",  resp = "y2", lb = 0),
		prior(normal(0.7, 0.7), nlpar = "A",  resp = "y2", lb = 0),
		prior_string(prior_r0, nlpar = "r0", resp = "y2", lb = 0)),
	chains = chains,
	iter = iter,
	cores = cores,
	control = list(
		adapt_delta = adapt_delta,
		max_treedepth = max_treedepth),
	backend = "cmdstanr",
	file = paste0("processed_data/Bayes_mich_men_",
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


### Bayesian parameterization following an exponential decay function

Bef_bayes_exp_decay<- function(data,
                                y_1,
                                y_2,
                                x,
                                group,
                                chains = 4,
                                iter = 4000,
                                cores = 4,
                                adapt_delta = 0.95,
                                max_treedepth = 12) {
	df <- data[, c(y_1, y_2, x, group)]
	names(df) <- c("y1", "y2", "x", "grp")

	df <- df[complete.cases(df[, c("x", "grp")]), ]
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

