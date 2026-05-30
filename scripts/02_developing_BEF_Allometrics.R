getwd()
setwd('/Users/hyli0001/wrd/b/Dynamic_allometrics/')
# system('ls -alt ../processed_data')


# fit_mm <- readRDS("processed_data/Bayes_mich_men_ft1.forest_type_1ch_20it_20260530_165427.rds")
# fit_exp <- readRDS("processed_data/Bayes_exp_decay_ft1.forest_type_1ch_20it_20260530_165427.rds")

# hist(as.data.frame(fit_mm)$b_y1_L_Intercept)
# plot(pred[,'Estimate','y1'], fit_mm$data$y1)
# plot(pred[,'Estimate','y2'], fit_mm$data$y1)

bp<-read.table('processed_data/plot_biomass.txt',sep='\t',head=TRUE)
head(bp)
##The number of spcies
sum(!is.na(unique(bp$sp_code)))
unique(bp$sp_code)
unique(bp$ft1.forest_type)
sum(!is.na(unique(bp$PFT)))

exp_decay<- function(x, L, A, k) {
  L + A * exp(-k * x)}

mich_men <- function(x,L,A,r0) {
  L + A / (1 + x/r0)}


par(mfrow=c(2,2))
for (i in unique(bp$ft1.forest_type)){
	df1<-bp[bp$ft1.forest_type==i,]

lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.shape,df1,ylim=c(1,lmx),main=i)
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1)
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1)
	points(beft.st~ sdi,df1,col=2)

	plot(befa.st~dist.mmd,ylim=c(1,lmx),df1)
	points(beft.st~ dist.mmd,df1,col=2)	
	
	}
unique(bp[bp$Family=='','sp_code'])
unique(bp[bp$sp_code=='LL','PFT'])
par(mfrow=c(2,2))
for (i in unique(bp$PFT)){
	df1<-bp[bp$PFT==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ cv,df1,ylim=c(1,lmx),main=i)
	points(beft.st~ cv,df1,col=2)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1)
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1)
	points(beft.st~ sdi,df1,col=2)
	
	plot(befa.st~dist.mmd,ylim=c(1,lmx),df1)
	points(beft.st~ dist.mmd,df1,col=2)	
	
	}

non.mod.mm<-nls(befa.st~mich_men(sdi,L,A,r0),bp[bp$Family=='Pinaceae',],start=c(L=0.3,A=1.2,r0=2.6))

non.mod.exp<-nls(befa.st~exp_decay(sdi,L,A,k),bp[bp$Family=='Pinaceae',],start=c(L=1.3,A=1.4,k=.02))


non.mod.mm<-nls(befa.st~mich_men(sdi,L,A,r0),bp[bp$PFT=='ENF',],start=c(L=0.3,A=1.2,r0=2.6))
non.mod.exp<-nls(befa.st~exp_decay(sdi,L,A,k),bp[bp$PFT=='ENF',],start=c(L=1.3,A=1.4,k=.02))


plot(befa.st~sdi,bp[bp$PFT=='ENF',],col=0)
for(i in 1:7){
	sp_c<-unique(bp[bp$PFT=='ENF','sp_code'])[i]
points(befa.st~sdi,bp[bp$sp_code==sp_c&bp$PFT=='ENF',],bg=i,pch=21)}
df00<-bp[bp$PFT=='ENF',]
sdi<-as.data.frame(seq(range(df00$sdi)[1],range(df00$sdi)[2],1))
colnames(sdi)<-'sdi'
bp$befa.a.pre<-predict(non.mod.mm,bp)
sdi$befa.a.pre<-predict(non.mod.mm,sdi)

summary(lm(befa.a.pre~befa.st,bp))
lines(befa.a.pre~sdi,sdi,lwd=1.1,col=4)

unique(bp[bp$PFT=='ENF','sp_code'])[1]
plot(befa.st~sdi,bp[bp$sp_code=='PK',])
points(befa.st~sdi,bp[bp$sp_code=='PK'&bp$stand_id==1005,],col=2)

plot(befa.st~sdi,bp[bp$Family=='Pinaceae',],col=0)
for(i in 1:6){
	sp_c<-unique(bp[bp$Family=='Pinaceae','sp_code'])[i]
points(befa.st~sdi,bp[bp$sp_code==sp_c&bp$Family=='Pinaceae',],bg=i,pch=21)}
points(befa.st~sdi,bp[bp$sp_code=='LL'&bp$Family=='Pinaceae',],bg='white',pch=22,col=4)
df00<-bp[bp$Family=='Pinaceae',]
sdi<-as.data.frame(seq(range(df00$sdi)[1],range(df00$sdi)[2],1))
colnames(sdi)<-'sdi'
bp$befa.a.pre<-predict(non.mod.mm,bp)
summary(lm(befa.a.pre~befa.st,bp))
lines(befa.a.pre~sdi,sdi,lwd=1.1,col=4)


par(mfrow=c(2,3))
for (i in unique(bp$Family)){
	df1<-bp[bp$Family==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1,main=i,xlim=c(0,.6))
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~ wb.shape,ylim=c(1,lmx),df1,main=i,xlim=c(0,12))
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1,main=i,xlim=c(0,700))
	points(beft.st~sdi,df1,col=2)	
	}


par(mfrow=c(2,3))
for (i in unique(bp$Genus)){
	df1<-bp[bp$Genus==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1,main=i,xlim=c(0,.6))
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~ wb.shape,ylim=c(1,lmx),df1,main=i,xlim=c(0,12))
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1,main=i,xlim=c(0,700))
	points(beft.st~ sdi,df1,col=2)	
	
	}

par(mfrow=c(2,3))
for (i in unique(bp$sp_code)[1:22]){
	df1<-bp[bp$sp_code==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1,main=i,xlim=c(0,.6))
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~ wb.shape,ylim=c(1,lmx),df1,main=i,xlim=c(0,12))
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1,main=i,xlim=c(0,700))
	points(beft.st~ sdi,df1,col=2)	
	}

plot(wb.scale~sdi,bp)











library(brms)
library(cmdstanr)

source('scripts/00_bayesian_functions.R')
