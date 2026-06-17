getwd()
setwd('/Users/hyli0001/wrd/b/Dynamic_allometrics/')
# system('ls -alt ../')
# bp<-read.table('processed_data/plot_biomass.txt',sep='\t',head=TRUE)
# baye.var<-c("stand_id","sp_code","plot","age","d","h","Bst","Bbr","Bf","Bcr","Vst","Vst.5","ba","befa.v","befa.st","beft.v","befr.st","beft.st","bwd","PFT","Genus","Family","ftp","ft1.forest_type","ft1.dominant_prop","sdi","sdi.1","sdi.2","sdi.3","sdi_max","rsd","rsd.1","rsd.2")
# write.table(bp[,baye.var],'processed_data/plot_biomass_bayes.txt',quote=FALSE,sep='\t',row.names=FALSE)
system('ls -alt bayes_outputs')
baydata<-read.table('processed_data/plot_biomass_bayes.txt',sep='\t',head=TRUE)
library(brms)

ftp_sp<-readRDS("bayes_outputs/exp_decay_ftp_sp_code_rsd_4chn_4000itr_4cor_0.99del_15depth.rds")
ftp<-readRDS("bayes_outputs/exp_decay_ftp_rsd_4chn_4000itr_4cor_0.99del_15depth.rds")
sp_c<-readRDS("bayes_outputs/exp_decay_sp_code_rsd_4chn_4000itr_4cor_0.99del_15depth.rds")
pft_sp<-readRDS("bayes_outputs/exp_decay_PFT_sp_code_rsd_4chn_4000itr_4cor_0.99del_15depth.rds")
pft<-readRDS("bayes_outputs/exp_decay_PFT_rsd_4chn_4000itr_4cor_0.99del_15depth.rds")


"exp_decay_log_partial_pool_ftp_sp_code_kdepth0_rsd_student_4chn_4000itr_4cor_0.99del_15depth.rds"
"exp_decay_log_partial_pool_ftp_sp_code_kdepth1_rsd_student_4chn_4000itr_4cor_0.99del_15depth.rds"
"exp_decay_log_partial_pool_ftp_sp_code_kdepth2_rsd_student_4chn_4000itr_4cor_0.99del_15depth.rds"

ftp.sp0<-readRDS("bayes_outputs/xp.ftp.sp_code.coef0.stud.rsd.4.4k.99.15.rds")
ftp.sp0

ftp.sp0n<-readRDS("bayes_outputs/exp_decay_log_partial_pool_ftp_sp_code_kdepth0_rsd_student_4chn_4000itr_4cor_0.99del_15depth.rds")
ftp.sp0n
source("scripts/00_bayesian_functions.R")
pp_check(ftp.sp0n,resp='y1')
pp_check(ftp.sp0,resp='z1')
bayes_backtransform_response("z1",yrep_z1)

ftp.sp0
summary(ftp_sp)
summary(pft_sp)
summary(ftp)
summary(pft)
summary(sp_c)

ranef(ft_sp)$h2[, , "y1_L_Intercept"]
coef(ft_sp)



h1_cf <- coef(ft_sp)$h1
h1_param <- data.frame(
  ft1_type = dimnames(h1_cf)[[1]],
  y1_L = h1_cf[, "Estimate", "y1_L_Intercept"],
  y1_A = h1_cf[, "Estimate", "y1_A_Intercept"],
  y1_k = h1_cf[, "Estimate", "y1_k_Intercept"],
  y2_L = h1_cf[, "Estimate", "y2_L_Intercept"],
  y2_A = h1_cf[, "Estimate", "y2_A_Intercept"],
  y2_k = h1_cf[, "Estimate", "y2_k_Intercept"]
)
h1_param$y_L <- h1_param$y1_L + h1_param$y2_L
h1_param$y_x0 <- (h1_param$y1_L + h1_param$y1_A) + (h1_param$y2_L + h1_param$y2_A)
h1_param$y_A <- h1_param$y1_A + h1_param$y2_A


h2_cf <- coef(ft_sp)$h2
h2_param <- data.frame(
  sp_code = dimnames(h2_cf)[[1]],
  y1_L = h2_cf[, "Estimate", "y1_L_Intercept"],
  y1_A = h2_cf[, "Estimate", "y1_A_Intercept"],
  y1_k = h2_cf[, "Estimate", "y1_k_Intercept"],
  y2_L = h2_cf[, "Estimate", "y2_L_Intercept"],
  y2_A = h2_cf[, "Estimate", "y2_A_Intercept"],
  y2_k = h2_cf[, "Estimate", "y2_k_Intercept"]
)
h2_param$y_L <- h2_param$y1_L+h2_param$y2_L
h2_param$y_x0 <- (h2_param$y1_L+h2_param$y1_A)+(h2_param$y2_L+h2_param$y2_A)
h2_param$y_A <- h2_param$y1_A + h2_param$y2_A
head(baydata)
range(baydata$d[baydata$PFT=='ENF'])
baydata[baydata$d<0.05,]



plot(befa.st~sdi,baydata[baydata$sp_code=='PDe',],ylim=c(1,6),xlim=c(0,500))
points(befa.st~sdi,baydata[baydata$sp_code=='PD',],col=4)
points(23,5.7,col=2)
points(23,4.69,col=2)


par(mfrow=c(1,3))
x<-'rsd'
y1<-'befa.st'
y2<-'beft.st'
for (i in unique(baydata$ft1.forest_type)){
	dft <- baydata[baydata$ft1.forest_type == i, ]

plot(dft[[x]],dft[[y2]], col = 0, xlab = "Relative stand density", ylab = "Biomass expansion factors",main = i,ylim = c(0.8,4),xlim=c(0,1.1))

	sp_list<-unique(dft$sp_code)
	for (j in seq_along(sp_list)){
		sp <- sp_list[j]
		dsub<-dft[dft$sp_code==sp,]
		points(dsub[[x]], dsub[[y1]], pch = 21, bg = j,cex=0.75, col=0)
		points(dsub[[x]], dsub[[y2]], pch = 2, col = j,lwd=0.5,cex=0.75)
		}

psub <- h1_param[h1_param$ft1_type == i, ]
xseq <- seq(min(dft[[x]], na.rm = TRUE),max(dft[[x]], na.rm = TRUE),length.out = 200)
ypre1 <- psub$y1_L + psub$y1_A * exp(-psub$y1_k * xseq)
ypre2 <- psub$y_L + psub$y_A * exp(-psub$y_x0 * xseq)

lines(xseq, ypre1, col = 1, lwd = 2)
lines(xseq, ypre2, col = 2, lwd = 2, lty = 2)
}

###################### PARAMETER DISTRIBUTION SHOULD BE MADE


############################## SPECIES TEST NEEDS TO BE DONE

for (i in seq_along(sp_list)) {
  sp <- sp_list[i]

  dsub <- dft[dft$sp_code == sp, ]
  psub <- h2_param[h2_param$sp_code == sp, ]

  if (nrow(psub) == 0) next

  points(dsub[[x]], dsub[[y1]], bg = i, pch = 21)
  points(dsub[[x]], dsub[[y2]], col = i, pch = 2)

  xseq <- seq(min(dsub[[x]], na.rm = TRUE),
              max(dsub[[x]], na.rm = TRUE),
              length.out = 200)

  ypre1 <- psub$y1_L + psub$y1_A * exp(-psub$y1_k * xseq)

  ypre2 <- psub$y2_L + psub$y2_A * exp(-psub$y2_k * xseq)

  lines(xseq, ypre1, col = i, lwd = 1.5)
  lines(xseq, ypre2, col = i, lwd = 1.5,lty=2)

}

legend("topright", legend = sp_list,
       col = seq_along(sp_list), pch = 21, lwd = 1.5)


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
	plot(befa.st~ rsd,df1,ylim=c(1,lmx),main=i)
	points(beft.st~rsd,df1,col=2)

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
