bayes_prior_normal <- function(mean, sd) {
	list(dist = "normal", mean = mean, sd = sd)
}

bayes_prior_exponential <- function(rate) {
	list(dist = "exponential", rate = rate)
}

bayes_shape_parameter <- function(pooling = c("full", "k_depth", "none"),
                                  sd_prior = bayes_prior_exponential(1)) {
	pooling <- match.arg(pooling)
	list(pooling = pooling, sd_prior = sd_prior)
}

bayes_model_shape <- function(name,
                              abbr = name,
                              mean,
                              parameters,
                              priors,
                              response_priors = list(sigma = bayes_prior_exponential(1))) {
	structure(
		list(
			name = name,
			abbr = abbr,
			mean = mean,
			parameters = parameters,
			priors = priors,
			response_priors = response_priors
		),
		class = "bayes_model_shape"
	)
}

bayes_model_shapes <- list(
	exp_decay = bayes_model_shape(
		name = "exp_decay",
		abbr = "xp",
		mean = "exp(logL) + exp(logA) * exp(-exp(logk) * x)",
		parameters = list(
			logL = bayes_shape_parameter("full"),
			logA = bayes_shape_parameter("full"),
			logk = bayes_shape_parameter("k_depth")
		),
		priors = list(
			z1 = list(
				logL = bayes_prior_normal(log(0.2), 0.8),
				logA = bayes_prior_normal(log(0.6), 0.7),
				logk = bayes_prior_normal(log(1.0), 0.8)
			),
			z2 = list(
				logL = bayes_prior_normal(log(0.2), 0.8),
				logA = bayes_prior_normal(log(0.7), 0.7),
				logk = bayes_prior_normal(log(1.0), 0.8)
			)
		)
	),

	michaelis_menten = bayes_model_shape(
		name = "michaelis_menten",
		abbr = "mm",
		mean = "exp(logL) + exp(logA) / (1 + x / exp(logr0))",
		parameters = list(
			logL = bayes_shape_parameter("full"),
			logA = bayes_shape_parameter("full"),
			logr0 = bayes_shape_parameter("k_depth")
		),
		priors = list(
			z1 = list(
				logL = bayes_prior_normal(log(0.2), 0.8),
				logA = bayes_prior_normal(log(0.6), 0.7),
				logr0 = bayes_prior_normal(log(1.0), 0.8)
			),
			z2 = list(
				logL = bayes_prior_normal(log(0.2), 0.8),
				logA = bayes_prior_normal(log(0.7), 0.7),
				logr0 = bayes_prior_normal(log(1.0), 0.8)
			)
		)
	),

	linear = bayes_model_shape(
		name = "linear",
		abbr = "ln",
		mean = "alpha + beta * x",
		parameters = list(
			alpha = bayes_shape_parameter("full"),
			beta = bayes_shape_parameter("k_depth")
		),
		priors = list(
			z1 = list(
				alpha = bayes_prior_normal(0, 1),
				beta = bayes_prior_normal(0, 1)
			),
			z2 = list(
				alpha = bayes_prior_normal(0, 1),
				beta = bayes_prior_normal(0, 1)
			)
		)
	)
)

bayes_model_shape_aliases <- c(
	exp = "exp_decay",
	exponential = "exp_decay",
	exponential_decay = "exp_decay",
	mm = "michaelis_menten",
	michaelis = "michaelis_menten",
	mich_menten = "michaelis_menten"
)

bayes_get_model_shape <- function(shape) {
	if (inherits(shape, "bayes_model_shape")) {
		return(shape)
	}

	name <- tolower(gsub("-", "_", trimws(as.character(shape)[1])))
	if (name %in% names(bayes_model_shape_aliases)) {
		name <- bayes_model_shape_aliases[[name]]
	}
	if (!name %in% names(bayes_model_shapes)) {
		stop(
			"Unknown model shape: ", shape,
			". Available shapes: ", paste(names(bayes_model_shapes), collapse = ", ")
		)
	}
	bayes_model_shapes[[name]]
}

bayes_active_model_shape <- "exp_decay"
