# renv::restore('..')
# renv::status('..')

source('scripts/00_basic_functions.R')

#################################
######	Loading dataset	 ########
#################################

raw.tb<-readRDS('raw_data/rawdata.rds')
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
ts<-merge(tsr,pft_ts,by='species_kr')
ts$PFTg<-paste(ts$PFT,ts$PFT_grow,sep='_')
pft_ts$PFTg<-paste(pft_ts$PFT,pft_ts$PFT_grow,sep='_')
setdiff(unique(ts$species),unique(pft_ts$species))

### By functional types
split_sp_pft<-split(ts,ts$stand_id)
divers.metrics<-lapply(split_sp_pft, function(df){
	multilevel_diversity(
		df$species,
		df$d,
		pft_ts,
		trait_group=c("PFT", "PFTg",'Family','Genus'))
		})

unlst <- lapply(divers.metrics, function(x) {
  c(
    sp.dom.p = unlist(x$dominant_prop),
    ft1 = unlist(x$Forest_type),
    ft2 = unlist(x$PFT_diversity),
    ft3 = unlist(x$PFTg_diversity),
    fml = unlist(x$Family_diversity),
    gns = unlist(x$Genus_diversity),
    sp = unlist(x$species_diversity),
    Gini_DBH = x$Gini_DBH,
    RaoQ_DBH = x$RaoQ_DBH
  )
})
dvst<-as.data.frame(do.call(rbind,unlst))
dvst$stand_id<-rownames(dvst)

### Top 3 species-level basal area and abundance per stand

sp_struct <- lapply(names(divers.metrics), function(id) {
  
  x <- divers.metrics[[id]]
  
  data.frame(
    stand_id = id,
    species  = names(x$species_basal_area),
    sp_ba    = as.numeric(x$species_basal_area),
    sp_n     = as.numeric(x$species_n),
    stringsAsFactors = FALSE
  )
})

sp_struct <- do.call(rbind, sp_struct)
sp_struct <- sp_struct[order(sp_struct$stand_id, -sp_struct$sp_ba), ]

sp_struct$rank <- ave(
  sp_struct$sp_ba,
  sp_struct$stand_id,
  FUN = function(z) rank(-z, ties.method = "first")
)

sp_top3 <- sp_struct[sp_struct$rank <= 3, ]

sp_top3_wide <- reshape(
  sp_top3[, c("stand_id", "rank", "species", "sp_ba", "sp_n")],
  idvar = "stand_id",
  timevar = "rank",
  direction = "wide"
)

wb.fit<-lapply(split(ts$d,ts$stand_id),weibull.para)
wb.para<-as.data.frame(do.call(rbind,wb.fit))
wb.para$stand_id<-rownames(wb.para)
rownames(wb.para)<-NULL

stand_sp<-unique(bh[,c('stand_id','sp_code')])
stand_sp<-stand_sp[!(stand_sp$stand_id==1001&stand_sp$sp_code=='CS'),]

ts.std<-merge(merge(merge(merge(merge(merge(plt,dvst,by='stand_id'),merge(collap(ts, d ~ stand_id + plot_size, FUN = list(m = function(x) sqrt(fmean(x^2)),med = function(x) sqrt(fmedian(x^2)),sd=fsd,n=fnobs)),collap(ts,h~stand_id+plot_size,FUN=list(m=fmean,med=fmedian,sd=fsd,n=fnobs)),by=c('stand_id','plot_size')),by='stand_id'),wb.para,by='stand_id'),sp_top3_wide,by='stand_id',all.x=TRUE),stand_sp,by='stand_id',all.x=TRUE),pft_h[, c("sp_code", "species")],by='sp_code',all.x=TRUE)

ts.std$std<-ts.std$n.d/ts.std$plot_size

ts.std$sp_ba.ha.1<-ts.std$sp_ba.1/ts.std$plot_size
ts.std$sp_std.1<-ts.std$sp_n.1/ts.std$plot_size
ts.std$sp_m.d.1<-sqrt((ts.std$sp_ba.ha.1/ts.std$sp_std.1)/pi)*2

ts.std$sp_ba.ha.2<-ts.std$sp_ba.2/ts.std$plot_size
ts.std$sp_std.2<-ts.std$sp_n.2/ts.std$plot_size
ts.std$sp_m.d.2<-sqrt((ts.std$sp_ba.ha.2/ts.std$sp_std.2)/pi)*2

ts.std$sp_ba.ha.3<-ts.std$sp_ba.3/ts.std$plot_size
ts.std$sp_std.3<-ts.std$sp_n.3/ts.std$plot_size
ts.std$sp_m.d.3<-sqrt((ts.std$sp_ba.ha.3/ts.std$sp_std.3)/pi)*2

ts.std$dist.mmd<-ts.std$med.d-ts.std$m.d
ts.std$sdi<-SDI(ts.std$std,ts.std$m.d*100,1.605)
ts.std$sdi.1<-SDI(ts.std$sp_std.1,ts.std$sp_m.d.1*100,1.605)
ts.std$sdi.2<-SDI(ts.std$sp_std.2,ts.std$sp_m.d.2*100,1.605)
ts.std$sdi.3<-SDI(ts.std$sp_std.3,ts.std$sp_m.d.3*100,1.605)

write.table(ts.std,'processed_data/stand_structure.txt',quote=FALSE,sep='\t',row.names=FALSE)


#################################
###### Analyzing NFI data ######
#################################
nfi<-read.table('original_files/KNFI/KNFI.txt',sep='\t',head=TRUE)

unique(nfi$sp.dom)[is.na(match(unique(nfi$sp.dom), pft_h$수종))]
pft_h$수종[is.na(match(pft_h$수종,unique(nfi$sp.dom)))]
nfi$sp.bef<-nfi$sp.dom
nfi$sp.bef[nfi$sp.dom=='편백']<-'편백나무'
nfi$sp.bef[nfi$sp.dom=='사시나무']<-'현사시나무'
nfi$sp.bef[nfi$sp.dom=='소나무']<-'중부지방소나무'
nfi$sp.bef[nfi$admin.1%in%c('강원도','강원특별자치도')&nfi$sp.dom=='소나무']<-'강원지방소나무'
pft_h$수종[is.na(match(pft_h$수종,unique(nfi$sp.bef)))]
nfi[nfi$sp.dom=='상수리나무',]
pft_h[pft_h$수종=='상수리',]


nfi_h<-merge(nfi[nfi$sp.bef%in%c(pft_h$수종),],pft_h,by.x='sp.bef',by.y='수종',all.x=TRUE)
nfi_h$d.qm<-sqrt((((nfi_h$ba.ha)/nfi_h$std.ha)/pi)*4)
nfi_h$sdi<-SDI(nfi_h$std.ha,nfi_h$d.qm*100,1.605)

head(nfi_h)

est_sdi_max<-function(df, group = "sp_code", sdi = "sdi", q = 0.99) {
  tapply(df[[sdi]], df[[group]], quantile, probs = q, na.rm = TRUE)
}
sdi_max_sp.nfi <- est_sdi_max(nfi_h, group = "sp_code", sdi = "sdi", q = 0.99)


#################################
### Analyzing harvested trees ###
#################################

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

bp<-merge(merge(bh[!is.na(bh$d),],pft_h[,c('species','sp_code','PFT','Genus','Family')],all.x=TRUE,by=c('sp_code')),ts.std,all.x=TRUE,by=c('stand_id','species','sp_code'))

bp$ba<-(bp$d/2)^2*pi
bp$befa.v<-(bp$Bst+bp$Bbr+bp$Bf)/bp$Vst
bp$befa.st<-(bp$Bst+bp$Bbr+bp$Bf)/bp$Bst
bp$beft.v<-(bp$Bst+bp$Bbr+bp$Bf+bp$Bcr)/bp$Vst
bp$beft.st<-(bp$Bst+bp$Bbr+bp$Bf+bp$Bcr)/bp$Bst
bp$befr.st<-bp$beft.st-bp$befa.st
bp$bwd<-bp$Bst/bp$Vst/1000
bp$d2h<-bp$d^2*bp$h
bp$cai<-(bp$Vst-bp$Vst.5)/5
bp$Bw<-(bp$Bst+bp$Bbr+bp$Bcr)
bp$Baw<-(bp$Bst+bp$Bbr)
bp$Ba<-(bp$Bst+bp$Bbr+bp$Bf)
bp$Bt<-(bp$Bst+bp$Bbr+bp$Bcr+bp$Bf)


sdi_max_sp.bp <- est_sdi_max(ts.std, group = "sp_code", sdi = "sdi", q = 0.99)
common_sp <- intersect(names(sdi_max_sp.bp), names(sdi_max_sp.nfi))
sdi_max_sp<-pmax(sdi_max_sp.bp[common_sp],sdi_max_sp.nfi[common_sp],na.rm=TRUE)

bp$sdi_max<-sdi_max_sp[as.character(bp$sp_code)]
bp$rsd<-bp$sdi/bp$sdi_max
bp$rsd.1<-bp$sdi.1/bp$sdi_max
bp$rsd.2<-bp$sdi.2/bp$sdi_max
bp$rsd.3<-bp$sdi.3/bp$sdi_max

bp$ftp<-bp$ft1.forest_type
bp$ftp[bp$ft1.forest_type=='mixed']<-'brd'
bp$ftp[bp$ft1.forest_type=='mono_B']<-'brd'
bp$ftp[bp$ft1.forest_type=='mono_N']<-'cnf'


write.table(bp,'processed_data/plot_biomass.txt',quote=FALSE,sep='\t',row.names=FALSE)




