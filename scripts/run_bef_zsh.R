bp<-read.table('processed_data/plot_biomass.txt',sep='\t',head=TRUE)

args <- commandArgs(trailingOnly = TRUE)

group_var <- args[1]
chains_n  <- as.numeric(args[2])
iter_n    <- as.numeric(args[3])
cores_n   <- as.numeric(args[4])
adapt_delta_n <- as.numeric(args[5])
max_treedepth_n   <- as.numeric(args[6])
x_var	<- args[7]

library(brms)
library(cmdstanr)

source("scripts/00_bayesian_functions.R")

fit_exp <- Bef_bayes_exp_decay(
  bp,
  y_1 = "befa.st",
  y_2 = "beft.st",
  x = x_var,
  group = group_var,
  chains = chains_n,
  iter = iter_n,
  cores = cores_n,
  adapt_delta = adapt_delta_n,
  max_treedepth = max_treedepth_n
)

fit_mm <- Bef_bayes_mich_men(
  bp,
  y_1 = "befa.st",
  y_2 = "beft.st",
  x = x_var,
  group = group_var,
  chains = chains_n,
  iter = iter_n,
  cores = cores_n,
  adapt_delta = adapt_delta_n,
  max_treedepth = max_treedepth_n
)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")


saveRDS(
  fit_mm,
  paste0(
    "processed_data/Bayes_mich_men_",
    group_var,
    "_",
    x_var,
    "_",
    chains_n, "chn_",
    iter_n, "itr_",
	cores_n, "cor_", 
    adapt_delta_n, "del_",
    max_treedepth_n, "depth",
    ".rds"
  )
)

saveRDS(
  fit_exp,
  paste0(
    "processed_data/Bayes_exp_decay_",
    group_var,
    "_",
    x_var,
    "_",
    chains_n, "chn_",
    iter_n, "itr_",
	cores_n, "cor_", 
    adapt_delta_n, "del_",
    max_treedepth_n, "depth",
    ".rds"
  )
)

tryCatch({
  ### fit_exp
  ### fit_mm
  ### saveRDS
  cat("\nSUCCESS\n")
  cat("Finished:", Sys.time(), "\n")
}, error = function(e) {
  cat("\nERROR OCCURRED\n")
  cat("Time:", Sys.time(), "\n")
  cat("Message:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})