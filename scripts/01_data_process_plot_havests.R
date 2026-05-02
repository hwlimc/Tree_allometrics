# renv::restore('..')
# renv::status('..')

getwd()
source('00_functions.R')
raw.tb<-readRDS('../raw_data/rawdata.rds')
str(raw.tb, max.level=2)

plt<-raw.tb$data$plot_info
pft_h<-raw.tb$data$pft_harvest
bh<-raw.tb$data$biomass_harvest
pft_ts<-raw.tb$data$pft_tree_survey
tsr<-raw.tb$data$tree_survey

setdiff(names(tsr),names(raw.tb$metadata$tree_survey$variables))
setdiff(unique(tsr$stand_id),unique(bh$stand_id))
setdiff(unique(tsr$species_kr),unique(pft_ts$species_kr))
setdiff(unique(bh$sp_code),unique(pft_h$sp_code))

#################################
### Analyzing stand structure ###
#################################
ts<-merge(tsr, pft_ts,by='species_kr')
ts$PFT3<-paste(ts$PFT1,ts$PFT2,sep='_')
pft_ts$PFT3<-paste(pft_ts$PFT1,pft_ts$PFT2,sep='_')

setdiff(unique(ts$species),unique(pft_ts$species))

### By functional types
split_sp_pft<-split(ts,ts$stand_id)
divers.metrics<-lapply(split_sp_pft, function(df){
	multilevel_diversity(
		df$species,
		df$d,
		pft_ts,
		trait_group=c("PFT1", "PFT3"))}
		)


unlst<-lapply(divers.metrics, function(x) {
  c(dom.sp.prop = unlist(x$dominant_prop),
	species = unlist(x$species_diversity),
	PFT1 = unlist(x$PFT1_diversity),
	PFT3 = unlist(x$PFT3_diversity))})
dvst<-as.data.frame(do.call(rbind,unlst))
dvst$stand_id<-rownames(dvst) 

wb.fit<-lapply(split(ts$d,ts$stand_id),weibull.para)
wb.para<-as.data.frame(do.call(rbind,wb.fit))
wb.para$stand_id<-rownames(wb.para)
rownames(wb.para)<-NULL

ts.std<-merge(merge(plt,dvst,by='stand_id'),merge(collap(ts,d+h~stand_id+plot_size,FUN=list(m=fmean,med=fmedian,sd=fsd,n=fnobs)),wb.para,by='stand_id'),by='stand_id')

ts.std$cv<-ts.std$sd.d/ts.std$m.d
ts.std$std<-ts.std$n.d/ts.std$plot_size
ts.std$dist.mmd<-ts.std$med.d-ts.std$m.d
ts.std$sdi<-SDI(ts.std$std,ts.std$m.d*100,1.605)
write.table(ts.std,'../processed_data/stand_structure.txt',quote=FALSE,sep='\t',row.names=FALSE)


pft_h<-raw.tb$data$pft_harvest
bh<-raw.tb$data$biomass_harvest

bh$ba<-(bh$d/2)^2*pi
bh$befa.v<-(bh$Bst+bh$Bbr+bh$Bf)/bh$Vst
bh$befa.st<-(bh$Bst+bh$Bbr+bh$Bf)/bh$Bst
bh$beft.v<-(bh$Bst+bh$Bbr+bh$Bf+bh$Bcr)/bh$Vst
bh$beft.st<-(bh$Bst+bh$Bbr+bh$Bf+bh$Bcr)/bh$Bst
bh$bwd<-bh$Bst/bh$Vst/1000

bh[bh$bwd>1,]
bh[bh$bwd<0.25,]
bh[bh$sp_code=='QS',]
qck<-bh[bh$stand_id%in%c(17,18,46:49,59,62,67,70),]
bh[bh$stand_id==18&bh$plot==2&bh$no==5,c('d')]<-NA
bh[bh$stand_id%in%46:49,c('d')]<-NA
bh[bh$stand_id==59&bh$plot==1&bh$no==3,c('d')]<-NA
bh[bh$stand_id==62&bh$plot==13&bh$no==4,c('d')]<-NA
bh[bh$stand_id==67&bh$plot==6&bh$no==5,c('d')]<-NA
bh[bh$stand_id==70,c('d')]<-NA
bh[bh$stand_id%in%1024:1025,c('d')]<-NA

bp<-merge(merge(bh[!is.na(bh$d),],pft_h[,c('species','sp_code','PFT','Genus','Family')],all.x=TRUE,by=c('sp_code')),ts.std,all.x=TRUE,by=c('stand_id'))

bp$ba<-(bp$d/2)^2*pi
bp$befa.v<-(bp$Bst+bp$Bbr+bp$Bf)/bp$Vst
bp$befa.st<-(bp$Bst+bp$Bbr+bp$Bf)/bp$Bst
bp$beft.v<-(bp$Bst+bp$Bbr+bp$Bf+bp$Bcr)/bp$Vst
bp$beft.st<-(bp$Bst+bp$Bbr+bp$Bf+bp$Bcr)/bp$Bst
bp$bwd<-bp$Bst/bp$Vst/1000


write.table(bp,'../processed_data/plot_biomass.txt',quote=FALSE,sep='\t',row.names=FALSE)




