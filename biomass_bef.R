# 붉가시	QAt
# 강원지방소나무	PDe
# 구실잣밤	CS
# 굴참	QV
# 낙엽송	LL
# 리기다	PR
# 백합	LT
# 삼나무	CJ
# 상수리	QAc
# 서어	CL
# 신갈	QM
# 자작	BP
# 잣	PK
# 졸참	QS
# 중부지방소나무	PD
# 편백	CO
# 해송	PT
# 현사시	PTl

source('functions.R')

bm<-read.table("tree_harv.txt",header=TRUE)
colnames(bm)<-c('ext_no','sp','plot','no','age','d','h','Bst','Bbr','Bf','Bcr','wd','Vst','Vst.5')
bm$ba<-(bm$d/2)^2*pi
bm$bef1<-(bm$Bst+bm$Bbr+bm$Bf)/bm$Bst
bm$bef2<-(bm$Bst+bm$Bbr+bm$Bf+bm$Bcr)/bm$Bst
bm$bef3<-(bm$Bst+bm$Bbr+bm$Bf+bm$Bcr)/bm$Vst

si1<-read.table("site_info_1.txt",header=TRUE)
si2<-read.table("site_info_2.txt",header=TRUE)
sti<-read.table("stand_info.txt",header=TRUE)
spk<-unique(sti[,c('sp.k','sp')])
spk[order(spk$sp.k),]

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
