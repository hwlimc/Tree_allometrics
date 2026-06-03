rm(list=ls())

library('collapse')
library('fitdistrplus')
SDI<-function(N.ha,D.cm,coef){
	((N.ha/2.4711)*(D.cm/25.4)^coef)}

diversity_metrics <- function(x) {
  p <- x / sum(x)
  p <- p[p > 0]
  
  H <- -sum(p * log(p))
  S <- length(p)
  
  list(
    richness = S,
    shannon = H,
    simpson = 1 - sum(p^2),
    evenness = ifelse(S > 1, H / log(S), NA)
  )
}
gini_numeric <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < 2 || mean(x, na.rm = TRUE) == 0) return(NA_real_)
  
  sum(abs(outer(x, x, "-"))) / (2 * n^2 * mean(x, na.rm = TRUE))
}


rao_q_dbh <- function(species, d, weight = "ba") {
  
  ba <- pi * (d / 2)^2
  
  sp_ba <- tapply(ba, species, sum, na.rm = TRUE)
  sp_d  <- tapply(d, species, mean, na.rm = TRUE)
  
  common_sp <- intersect(names(sp_ba), names(sp_d))
  sp_ba <- sp_ba[common_sp]
  sp_d  <- sp_d[common_sp]
  
  if (length(sp_d) < 2) return(0)
  
  if (weight == "ba") {
    p <- sp_ba / sum(sp_ba, na.rm = TRUE)
  } else {
    sp_n <- table(species)[common_sp]
    p <- sp_n / sum(sp_n, na.rm = TRUE)
  }
  
  dist_mat <- abs(outer(sp_d, sp_d, "-"))
  
  sum(outer(p, p) * dist_mat, na.rm = TRUE)
}


multilevel_diversity <- function(species, d, trait_table, trait_group) {
  
  ba <- pi * (d / 2)^2
  
  species_ba <- tapply(ba, species, sum, na.rm = TRUE)
  species_n  <- table(species)
  
  output <- list(
    species_basal_area = species_ba,
    species_n = species_n,
    dom.sp.prop = max(species_ba, na.rm = TRUE) / sum(species_ba, na.rm = TRUE),
    species_diversity = diversity_metrics(species_ba),
    Gini_DBH = gini_numeric(d),
    RaoQ_DBH = rao_q_dbh(species, d, weight = "ba")
  )
  
  for (coln in trait_group) {
    
    group <- trait_table[[coln]][match(species, trait_table$species)]
    group_ba <- tapply(ba, group, sum, na.rm = TRUE)
    
    output[[paste0(coln, "_basal_area")]] <- group_ba
    output[[paste0(coln, "_diversity")]] <- diversity_metrics(group_ba)
  }
  
  prop <- tapply(
    ba,
    substr(trait_table$PFT[match(species, trait_table$species)], 2, 2),
    sum
  )
  
  prop <- prop / sum(prop)
  
  B <- prop["B"]
  N <- prop["N"]
  
  B[is.na(B)] <- 0
  N[is.na(N)] <- 0
 ### Forest type as needle leaves, broad leaves, and mixture between the two 
  output[["Forest_type"]] <- list(
    forest_type = c("mixed", "mono_B", "mono_N")[1 + (B >= .75) + 2 * (N >= .75)],
    dominant_prop = max(B, N)
  )
  
  output
}


weibull.para<-function(x){
	x<-x[!is.na(x)]
	fit<-fitdist(x,'weibull')
	c(wb=fit$estimate[1],
	wb=fit$estimate[2],
	se=fit$sd[1],
	se=fit$sd[2])}
