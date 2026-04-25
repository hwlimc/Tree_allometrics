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



write.table(ts.std,'../processed_data/stand_structure.txt',quote=FALSE,sep='\t')




bh$ba<-(bh$d/2)^2*pi
bh$befa.v<-(bh$Bst+bh$Bbr+bh$Bf)/bh$Vst
bh$befa.st<-(bh$Bst+bh$Bbr+bh$Bf)/bh$Bst
bh$beft.v<-(bh$Bst+bh$Bbr+bh$Bf+bh$Bcr)/bh$Vst
bh$beft.st<-(bh$Bst+bh$Bbr+bh$Bf+bh$Bcr)/bh$Bst

plot(Bst~Vst,bh)
abline(a=0,b=1000)
plot(wd~d,bh)
plot(Bbr~d,bh)
plot(Bst~d,bh)
plot(Vst~d,bh)


collap(ts,d+h~stand_id+plot_size,FUN=list(m=fmean,med=fmedian,sd=fsd,n=fnobs))

head(ts.std)



si1$ba<-(si1$dbh/200)^2*pi
sim1<-summaryBy(ba+h~ext_no+plot_size,si1,FUN=md)

sim2<-summaryBy(ba~ext_no+plot_size,si1,FUN=n)

sim2$std<-sim2$ba.n/sim2$plot_size
sim2$sp.kr<-sim2$sp
sim2$sp<-NULL
sim<-merge(merge(sim1,si2,by='ext_no'),sim2[,c('ext_no','std')],by='ext_no')
sim$rsd<-((sim$std/2.4711)*(sqrt(sim$ba.m/pi)*2*100/25.4)^1.605)/400

bms<-merge(bm,sim,by=c('ext_no'))
unique(bms$sp)


length(unique(bms$sp))

head(bms)

par(mfrow=c(3,6))
for(i in 1:18){
df<-bms[bms$sp==unique(bms$sp)[i],]
plot(bef2~rsd,df,bg=2,pch=21,ylim=c(1,max(df$bef2,na.rm=TRUE)),main=unique(bms$sp)[i])
points(bef1~rsd,df)}

plot(bef2~rsd,bms[bms$sp=='QAt',],bg=2,pch=21,ylim=c(0,4))
points(bef1~rsd,bms[bms$sp=='QAt',])
summary(lm(bef2~rsd,bms[bms$sp=='QAt',]))
head(bms)

plot(res1~rsd,bms[bms$sp=='QAt',],col=plot)
head(bms[bms$sp=='QAt',])
bms$res1<-NA
bms$pre<-NA
bms[bms$sp=='QAt',]$res<-residuals(lm(log(Bbr)~log(d)+rsd,bms[bms$sp=='QAt',]))

bms[bms$sp=='QAt',]$pre<-exp(predict(lm(log(Bbr)~log(d),bms[bms$sp=='QAt',]),bms[bms$sp=='QAt',]))
bms[bms$sp=='QAt',]$res1<-bms[bms$sp=='QAt',]$Bbr-bms[bms$sp=='QAt',]$pre


plot(bef2~age,bms,ylim=c(1,2.35),bg=2,pch=21)
points(bef1~age,bms,ylim=c(1,2.35))
