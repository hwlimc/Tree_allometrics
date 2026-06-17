getwd()
setwd('/Users/hyli0001/wrd/b/Dynamic_allometrics/')
# system('ls -alt ../processed_data')
bp<-read.table('processed_data/plot_biomass.txt',sep='\t',head=TRUE)

## The number of spcies
sum(!is.na(unique(bp$sp_code)))
unique(bp$sp_code)
unique(bp$ftp)
unique(bp$PFT)

exp_decay<- function(x, L, A, k) {
  L + A * exp(-k * x)
}

mod.eb<-nls(befa.st~exp_decay(rsd,L,A,k),bp[bp$PFT=='EBF',],start=c(L=2,A=4,k=1))
mod.db<-nls(befa.st~exp_decay(rsd,L,A,k),bp[bp$PFT=='DBF',],start=c(L=1.3,A=1.4,k=2))
mod.en<-nls(befa.st~exp_decay(rsd,L,A,k),bp[bp$PFT=='ENF',],start=c(L=1.3,A=1.4,k=2))
mod.dn<-nls(befa.st~exp_decay(rsd,L,A,k),bp[bp$PFT=='DNF',],start=c(L=1.2,A=4,k=1))

plot(befa.st~rsd,bp[bp$PFT=='EBF',])
unique(bp[bp$PFT=='EBF','sp_code'])


plot(befa.st~rsd,bp[bp$PFT=='EBF',])
points(befa.st~rsd,bp[bp$PFT=='EBF'&bp$sp_code==unique(bp[bp$PFT=='EBF','sp_code'])[1],],bg=1,pch=21)
points(befa.st~rsd,bp[bp$PFT=='EBF'&bp$sp_code==unique(bp[bp$PFT=='EBF','sp_code'])[7],],bg=7,pch=21)


unique(bp[bp$PFT=='EBF','sp_code'])[1] ## OK
unique(bp[bp$PFT=='EBF','sp_code'])[2] ## OK
unique(bp[bp$PFT=='EBF','sp_code'])[3] ## Not enough samples
unique(bp[bp$PFT=='EBF','sp_code'])[4] ## Opposite trend
unique(bp[bp$PFT=='EBF','sp_code'])[5] ## Not enough samples
unique(bp[bp$PFT=='EBF','sp_code'])[6] ## Not enough samples

unique(bp[bp$PFT=='DBF','sp_code'])

plot(befa.st~rsd,bp[bp$PFT=='DBF',])
points(befa.st~rsd,bp[bp$PFT=='DBF'&bp$sp_code==unique(bp[bp$PFT=='DBF','sp_code'])[9],],bg=4,pch=21)

unique(bp[bp$PFT=='DBF','sp_code'])[1] ## OK 
unique(bp[bp$PFT=='DBF','sp_code'])[2] ## OK
unique(bp[bp$PFT=='DBF','sp_code'])[3] ## OK
unique(bp[bp$PFT=='DBF','sp_code'])[4] ## Not enough samples
unique(bp[bp$PFT=='DBF','sp_code'])[5] ## Not enough samples
unique(bp[bp$PFT=='DBF','sp_code'])[6] ## No pattern
unique(bp[bp$PFT=='DBF','sp_code'])[7] ## OK
unique(bp[bp$PFT=='DBF','sp_code'])[8] ## OK



unique(bp[bp$PFT=='ENF','sp_code'])

plot(befa.st~rsd,bp[bp$PFT=='ENF',])
points(befa.st~rsd,bp[bp$PFT=='ENF'&bp$sp_code==unique(bp[bp$PFT=='ENF','sp_code'])[5],],bg=5,pch=21)

unique(bp[bp$PFT=='ENF','sp_code'])[1] ## OK
unique(bp[bp$PFT=='ENF','sp_code'])[2] ## OK
unique(bp[bp$PFT=='ENF','sp_code'])[3] ## OK
unique(bp[bp$PFT=='ENF','sp_code'])[4] ## OK
unique(bp[bp$PFT=='ENF','sp_code'])[5] ## OK
unique(bp[bp$PFT=='ENF','sp_code'])[6] ## OK
unique(bp[bp$PFT=='ENF','sp_code'])[7] ## OK




nrow(bp[bp$PFT=='EBF',])
summary(mod.enf)

mod.<-nls(befa.st~exp_decay(rsd,L,A,k),bp[bp$ft1.forest_type=='mono_B',],start=c(L=1.3,A=1.4,k=2))

anova(non.mod.exp1, non.mod.exp)
head(bp)


summary(non.mod.exp)

non.mod.exp1<-nls(befa.st~exp_decay(rsd,L,A,k),baydata[baydata$ft1.forest_type=='mono_B',],start=c(L=1.3,A=1.4,k=.2))
summary(non.mod.exp1)
anova(non.mod.exp1, non.mod.exp)



plot(befa.st~sdi,bp[bp$PFT=='ENF',],col=0)
for(i in 1:6){
	sp_c<-unique(bp[bp$PFT=='ENF','sp_code'])[i]
points(befa.st~sdi,bp[bp$sp_code==sp_c&bp$PFT=='ENF',],bg=i,pch=21)}
df00<-bp[bp$PFT=='ENF',]
sdi<-as.data.frame(seq(range(df00$sdi)[1],range(df00$sdi)[2],1))
colnames(sdi)<-'sdi'
bp$befa.a.pre<-predict(non.mod.mm,bp)
sdi$befa.a.pre<-predict(non.mod.mm,sdi)

summary(lm(befa.a.pre~befa.st,bp))
lines(befa.a.pre~sdi,sdi,lwd=1.1,col=4)




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
