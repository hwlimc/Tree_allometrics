bp<-read.table('processed_data/plot_biomass.txt',sep='\t',head=TRUE)
bp$d2h<-bp$d^2*bp$h
bp$cai<-(bp$Vst-bp$Vst.5)/5
bp$Bw<-(bp$Bst+bp$Bbr+bp$Bcr)
bp$Baw<-(bp$Bst+bp$Bbr)
bp$Ba<-(bp$Bst+bp$Bbr+bp$Bf)
bp$Bt<-(bp$Bst+bp$Bbr+bp$Bcr+bp$Bf)

args <- commandArgs(trailingOnly = TRUE)

group_var <- args[1]
chains_n  <- as.numeric(args[2])
iter_n    <- as.numeric(args[3])
cores_n   <- as.numeric(args[4])
adapt_delta_n <- as.numeric(args[5])
treedepth_n   <- as.numeric(args[6])

library(brms)
library(cmdstanr)

source("scripts/00_bayesian_functions.R")

fit_exp <- Bef_bayes_exp_decay(
  bp,
  y_1 = "befa.st",
  y_2 = "beft.st",
  x = "sdi",
  group = group_var,
  chains = chains_n,
  iter = iter_n,
  cores = cores_n,
  adapt_delta = adapt_delta_n,
  max_treedepth = treedepth_n
)

fit_mm <- Bef_bayes_mich_men(
  bp,
  y_1 = "befa.st",
  y_2 = "beft.st",
  x = "sdi",
  group = group_var,
  chains = chains_n,
  iter = iter_n,
  cores = cores_n,
  adapt_delta = adapt_delta_n,
  max_treedepth = treedepth_n
)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")


saveRDS(
  fit_exp,
  paste0(
    "processed_data/Bayes_exp_decay_",
    group_var,
    "_",
    chains_n, "ch_",
    iter_n, "it_",
    timestamp,
    ".rds"
  )
)

saveRDS(
  fit_mm,
  paste0(
    "processed_data/Bayes_mich_men_",
    group_var,
    "_",
    chains_n, "ch_",
    iter_n, "it_",
    timestamp,
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