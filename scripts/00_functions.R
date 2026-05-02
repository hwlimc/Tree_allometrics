rm(list=ls())
library('renv')
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

multilevel_diversity <- function(species, d, trait_table, trait_group)
{
	ba <- pi * (d/2)^2	# Species-level basal area
	species_ba <- tapply(ba, species, sum,na.rm=TRUE)

output <- list(
	species_basal_area = species_ba,
	dom.sp.prop = max(species_ba,na.rm=TRUE)/sum(species_ba,na.rm=TRUE),
	species_diversity = diversity_metrics(species_ba)
  )
  
  	for(coln in trait_group) {
		group <- trait_table[[coln]][match(species, trait_table$species)]
		group_ba <- tapply(ba, group, sum,na.rm=TRUE)
    
		output[[paste0(coln, "_basal_area")]] <- group_ba
		output[[paste0(coln, "_diversity")]] <- diversity_metrics(group_ba)
		}

	prop<-tapply(ba, substr(trait_table$PFT[match(species, trait_table$species)],2,2), sum)
	prop<-prop/sum(prop)
	
	B <- prop["B"]
	N <- prop["N"]
	B[is.na(B)] <- 0
	N[is.na(N)] <- 0
	
	
	output[['Forest_type']] <- list(
		forest_type = c("mixed","mono_B","mono_N")[1 + (B >= .75) + 2*(N >= .75)],
		dominant_prop = max(B,N))
output
}


weibull.para<-function(x){
	x<-x[!is.na(x)]
	fit<-fitdist(x,'weibull')
	c(wb=fit$estimate[1],
	wb=fit$estimate[2],
	se=fit$sd[1],
	se=fit$sd[2])}
