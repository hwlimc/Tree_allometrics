getwd()
setwd('/Users/hyli0001/wrd/b/Dynamic_allometrics/')
# system('ls -alt ../processed_data')
bp<-read.table('processed_data/plot_biomass.txt',sep='\t',head=TRUE)

# baye.var<-c("stand_id","sp_code","plot","age","d","h","Bst","Bbr","Bf","Bcr","Vst","Vst.5","ba","befa.v","befa.st","beft.v","beft.st","bwd","PFT","Genus","Family","ft1.forest_type","ft1.dominant_prop","sdi","sdi.1","sdi.2","sdi.3","sdi_max","rsd","rsd.1","rsd.2")
# write.table(bp[,baye.var],'processed_data/plot_biomass_bayes.txt',quote=FALSE,sep='\t',row.names=FALSE)


library(brms)


fit_mm<-readRDS("processed_data/Bayes_mich_men_ft1.forest_type_rsd.1_4chn_8000itr_4cor_0.99del_15depth.rds")
summary(fit_mm)
names(fit_mm$data)
exp_decay<- function(x, L, A, k) {
  L + A * exp(-k * x)}

mich_men <- function(x,L,A,r0) {
  L + A / (1 + x/r0)}


fit_exp<-readRDS("processed_data/Bayes_exp_decay_sp_code_rsd.1_4chn_8000itr_4cor_0.99del_15depth.rds")
summary(fit_exp)
pred<-fitted(fit_exp)
posterior(fit_exp)
names(fit_exp$data)
hist(as.data.frame(fit_exp)$b_y1_L_Intercept)
plot(pred[,'Estimate','y1'], fit_exp$data$y1,ylim=c(1,5),xlim=c(1,5))
plot(pred[,'Estimate','y2'], fit_exp$data$y2,ylim=c(1,5),xlim=c(1,5))

bp<-read.table('processed_data/plot_biomass.txt',sep='\t',head=TRUE)
head(bp)

source('scripts/00_bayesian_exp_decay.R')

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(brms)
library(cmdstanr)
fit <- Bef_bayes_exp_decay(
  data = bp,
  y_1 = "befa.st",
  y_2 = "beft.st",
  x = "rsd",
  hierarchy = c("ft1.forest_type","sp_code"),
  chains = 2,
  iter = 1000,
  cores = 2,
  adapt_delta = 0.95,
  max_treedepth = 12
)

plot(fit)


bppd<-bp[bp$sp_code%in%c('PDe','PD'),]

<-bp$sp_code[bp$sp_code=='PDe']

plot(befa.st~age, bppd)
points(beft.st~age, bppd[bppd$sp_code=='PDe',],col=2)
non.mod.exp<-nls(befa.st~exp_decay(rsd,L,A,k),bppd,start=c(L=11,A=12,k=3))

curve(exp_decay(x,summary(non.mod.exp)$coe[1],summary(non.mod.exp)$coe[2],summary(non.mod.exp)$coe[3]),add=TRUE)


bppd$befa.st.pre<-predict(non.mod.exp,bppd)
bppd$befa.st.res<-bppd$befa.st-bppd$befa.st.pre
curve(exp_decay(x,summary(non.mod.exp)$coe[1],summary(non.mod.exp)$coe[2],summary(non.mod.exp)$coe[3]),add=TRUE)
lines(befa.st.pre~rsd,rsd,lwd=1.1,col=4)
summary(lm(befa.st~befa.st.pre*sp_code,bppd))
plot(befa.st.res~befa.st.pre,bppd)
points(befa.st.res~befa.st.pre,bppd[bppd$sp_code=='PDe',],bg=2,pch=21)


unique(bp[,c('species','sp_code')])
head(bp)
##The number of spcies
sum(!is.na(unique(bp$sp_code)))
unique(bp$sp_code)
unique(bp$ft1.forest_type)
sum(!is.na(unique(bp$PFT)))



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

	plot(befa.st~rsd,ylim=c(1,lmx),df1)
	points(beft.st~ rsd,df1,col=2)	
	}
	

par(mfrow=c(2,2))
for (i in unique(bp$PFT)){
	df1<-bp[bp$PFT==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.shape,df1,ylim=c(1,lmx),main=i)
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1)
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1)
	points(beft.st~ sdi,df1,col=2)
	
	plot(befa.st~rsd,ylim=c(1,lmx),df1)
	points(beft.st~ rsd,df1,col=2)		
	}

par(mfrow=c(2,2))
for (i in unique(bp$sp_code)){
	df1<-bp[bp$sp_code==i,]
lmx<-max(c(df1$beft.st,df1$befa.st),na.rm=TRUE)
	plot(befa.st~ wb.shape,df1,ylim=c(1,lmx),main=i)
	points(beft.st~ wb.shape,df1,col=2)

	plot(befa.st~ wb.scale,ylim=c(1,lmx),df1)
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,ylim=c(1,lmx),df1)
	points(beft.st~ sdi,df1,col=2)
	
	plot(befa.st~rsd,ylim=c(1,lmx),df1)
	points(beft.st~ rsd,df1,col=2)		
	}



non.mod.exp<-nls(befa.st~exp_decay(rsd,L,A,k),bp[bp$PFT=='ENF',],start=c(L=1.3,A=1.4,k=.2))

plot(befa.st~rsd,bp[bp$PFT=='ENF',],col=0)
for(i in 1:7){
	sp_c<-unique(bp[bp$PFT=='ENF','sp_code'])[i]
points(befa.st~rsd,bp[bp$sp_code==sp_c&bp$PFT=='ENF',],bg=i,pch=21)}
df00<-bp[bp$PFT=='ENF',]
rsd<-as.data.frame(seq(range(df00$rsd)[1],range(df00$rsd)[2],0.001))
colnames(rsd)<-'rsd'
bp$befa.a.pre<-predict(non.mod.exp,bp)
rsd$befa.a.pre<-predict(non.mod.exp,rsd)
summary(lm(befa.a.pre~befa.st,bp))
lines(befa.a.pre~rsd,rsd,lwd=1.1,col=4)

non.mod.exp<-nls(befa.st~exp_decay(sdi,L,A,k),bp[bp$PFT=='ENF',],start=c(L=1.3,A=1.4,k=0.1))

plot(befa.st~sdi,bp[bp$PFT=='ENF',],col=0)
for(i in 1:7){
	sp_c<-unique(bp[bp$PFT=='ENF','sp_code'])[i]
points(befa.st~sdi,bp[bp$sp_code==sp_c&bp$PFT=='ENF',],bg=i,pch=21)}
df00<-bp[bp$PFT=='ENF',]
sdi<-as.data.frame(seq(range(df00$sdi)[1],range(df00$sdi)[2],0.001))
colnames(sdi)<-'sdi'
bp$befa.a.pre<-predict(non.mod.exp,bp)
sdi$befa.a.pre<-predict(non.mod.exp,sdi)
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
