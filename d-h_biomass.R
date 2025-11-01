source('~/Library/CloudStorage/OneDrive-Sverigeslantbruksuniversitet/wd/Rwork/functions.R')

bm<-read.table("tree_harv.txt",header=TRUE)
si1<-read.table("site_info_1.txt",header=TRUE)
si2<-read.table("site_info_2.txt",header=TRUE)
head(si1)
########### grouping slopes and find the effect of slope on H~D relationships (Fayolle et al. 2016; Pual et al. 2018)
head(si1,15)
si1$dbh<-si1$dbh/100
head(bm,15)

std<-aggregate(dbh~ext_no,data=si1,length)
std$plot.size<-aggregate(plot_size~ext_no,data=si1,mean)[,2]
std$qm.dbh<-aggregate(dbh~ext_no,data=si1,rmse)[,2]
std$sd<-std$dbh/std$plot.size
std<-unique(join(std,bm[,1:2],by=("ext_no"),type="inner"))
std<-unique(join(std,si2,by=("ext_no"),type="inner"))
colnames(std)<-c("ext_no","num_t","p_size","qm.dbh","sd","sp",colnames(si2)[2:6])
head(std)

length(std[std$sp%in%c('PD','PDe'),1])
plot(getMap(resolution="high"),xlim=c(125,130),ylim=c(33.2,38.5))
points(lati.~long.,std[std$sp%in%c('PD'),],pch=21,bg=1,cex=.5)
points(lati.~long.,std[std$sp%in%c('PDe'),],pch=21,bg=2,cex=.5)


m.age<-aggregate(age~ext_no,data=bm,mean)
colnames(m.age)[2]<-"m.age"
cv.age<-aggregate(age~ext_no,data=bm,cv)
colnames(cv.age)[2]<-"cv.age"
maxh1<-aggregate(h~ext_no,data=si1,max)
colnames(maxh1)[2]<-"maxh1"
maxh2<-aggregate(h~ext_no,data=bm,max)
colnames(maxh2)[2]<-"maxh2"
sp.k<-unique(si1[,1:2])
colnames(sp.k)[2]<-"sp.k"

s.info<-join(std,m.age,by=("ext_no"))
s.info<-join(s.info,cv.age,by=("ext_no"))
s.info<-join(s.info,maxh1,by=("ext_no"))
s.info<-join(s.info,maxh2,by=("ext_no"))
s.info<-join(s.info,sp.k,by=("ext_no"))
write.table(s.info,"stand_info.txt",sep="\t",quote=FALSE,row.names=FALSE)
n.stand<-aggregate(m.age~sp.k,data=s.info,length)


bm<-join(bm,std,by=("ext_no"))
bm$plot<-as.factor(bm$plot)
bm[,19]<-NULL
bm$Mab<-bm$Mst+bm$br+bm$leaf
bm$br<-ifelse(bm$br==0,bm$br+0.00001,bm$br)
bm$Mw<-bm$Mst+bm$br
bm$Mt<-bm$Mab+bm$cr
bm$rab<-bm$Mab/bm$Mt
bm$rbw<-bm$cr/bm$Mt
bm$abr<-bm$Mab/bm$cr
bm$BEFab<-bm$Mab/bm$Vst
bm$BEF<-bm$Mt/bm$Vst
bm$dV<-(bm$Vst-bm$Vst.5)/5
bm$n.sp<-as.factor(as.numeric(bm$sp))
bm$ln_dbh<-log(bm$dbh)
bm$ln_h<-log(bm$h)
bm$ln_Mst<-log(bm$Mst)
bm$ln_Mw<-log(bm$Mw)
bm$ln_br<-log(bm$br)
bm$ln_cr<-log(bm$cr)
bm$ln_leaf<-log(bm$leaf)
bm$ln_qm.dbh<-log(bm$qm.dbh)
bm$ln_sd<-log(bm$sd)
bm$d2h<-bm$dbh^2*bm$h
bm$ln_d2h<-log(bm$d2h)
bm$r.ba<-(bm$dbh^2)/(bm$qm.dbh^2)
head(bm,15)
species<-unique(bm$sp)
par(mfrow=c(3,8))
for (i in species){
	plot(dV~leaf,data=bm[bm$sp==i,])
	}

for (i in species){
	print(i)
	print(length(unique(bm[bm$sp==i,]$plot)))
	}

plot(Mab~cr,data=bm[bm$sp=="PD",])
plot(rab~leaf,data=bm[bm$sp=="PD",])
plot(rbw~leaf,data=bm[bm$sp=="PD",])
summary(lm(dV~leaf*lati.+alti.*leaf,data=bm[bm$sp=="PD",]))

# plot(qm.dbh~sd,data=bm)
# plot(ln_qm.dbh~ln_sd,data=bm)
# abline(lm(ln_qm.dbh~ln_sd,data=bm))

# most<-by(bm, bm$sp, function(bm) summary(lm(ln_Mst~ln_dbh,data=bm)))
# mobr<-by(bm, bm$sp, function(bm) summary(lm(ln_br~ln_dbh,data=bm)))
# mocr<-by(bm, bm$sp, function(bm) summary(lm(ln_cr~ln_dbh,data=bm)))

pine<-subset(bm,bm$sp=="PD")
h_model.1<-lm(ln_h~ln_dbh,data=pine)
h_model.2<-lm(ln_h~ln_dbh+plot,data=pine)
h_model.3<-lm(ln_h~ln_dbh+qm.dbh,data=pine)
h_model.4<-lm(ln_h~ln_dbh+r.ba,data=pine)
h_model.5<-lm(ln_h~ln_dbh+std,data=pine)
h_model.6<-lm(ln_h~ln_dbh+alti.,data=pine)
h_model.7<-lm(ln_h~ln_dbh+lati.,data=pine)
h_model.8<-lm(ln_h~ln_dbh+slope,data=pine)

s1<-summary(h_model.1)
s2<-summary(h_model.2)
s3<-summary(h_model.3)
s4<-summary(h_model.4)
s5<-summary(h_model.5)
s6<-summary(h_model.6)
s7<-summary(h_model.7)
s8<-summary(h_model.8)


anova(h_model.1)
anova(h_model.2)
anova(h_model.3)
anova(h_model.4)
anova(h_model.5)
anova(h_model.6)
anova(h_model.7)
anova(h_model.8)


rmdhp<-lme(ln_h~ln_dbh,random=~1|plot,data=pine)
summary(rmdhp)$
ranef(rmdhp)
anova(rmdhp,h_model.1)

gm<-pine$h-exp(predict(h_model.1))
sqgm<-sum(gm^2)
MSE1<-sqgm/length(gm)
rgm<-pine$h-exp(predict(rmdhp))
sqrgm<-sum(rgm^2)
MSE2<-sqrgm/length(rgm)
cbind(MSE1,MSE2)

plot(NA,xlim=c(0,0.36),ylim=c(2,21),xaxt="n", yaxt="n",ylab="",xlab="")
points(h~dbh,data=pine,pch=21,bg=1,cex=0.5)
for(i in 1:max(as.numeric(pine$plot))){
curve(exp(summary(rmdhp)$coef$fixed[1]+summary(rmdhp)$coef$random$plot[i])*(x^summary(rmdhp)$coef$fixed[2]),add=TRUE,col=i,xlim=c(min(pine[pine$plot==i,]$dbh),max(pine[pine$plot==i,]$dbh)))
}
curve(exp(summary(rmdhp)$coef$fixed[1])*(x^summary(rmdhp)$coef$fixed[2]),add=TRUE,col=1,lwd=2,xlim=c(min(pine$dbh),max(pine$dbh)))

axis (1, seq(0,0.4,by=0.1),tck=.02,mgp=c(0,.1,0))
axis (2, seq(0,22,by=5),tck=.02,mgp=c(0,.1,0))
title(xlab=expression(paste("Diameter at 1.2 m (m)")), ylab=expression(paste("Height (m)")), mgp=c(1.25,0,0))
mtext("Pinus densiflora",side=3,adj=0,padj=0.1,font=4)
mtext("South Korea",side=3,adj=1,padj=0.2,font=2)  


plot(h~dbh,data=pine)
A<-cbind(pine$dbh,exp(predict(rmdhp)))
A<-A[order(A[,1]),]
lines(A[,1],A[,2])

mo_dbh<-by(pine, pine$plot, function(pine) summary(lm(ln_Mst~ln_dbh,data=pine)))


coeff_d <- dlply(pine, .(plot), lm, formula=ln_Mst~ln_dbh)
cfd<-ldply(coeff_d, extractfun)
coeff_dh <- dlply(pine, .(plot), lm, formula=ln_Mst~ln_d2h)
cfdh<-ldply(coeff_dh, extractfun)

par(mfrow=c(1,2))
plot(ln_Mst~ln_dbh,data=pine,col=p.col)
for(i in 1:max(as.numeric(pine$plot))){
	curve(cfd[3][cfd$plot==i,]*x+cfd[2][cfd$plot==i,],col=i,xlim=c(min(pine[pine$plot==i,]$ln_dbh),max(pine[pine$plot==i,]$ln_dbh)),add=TRUE)
	}

plot(ln_Mst~ln_d2h,data=pine,col=p.col)
for(i in 1:max(as.numeric(pine$plot))){
	curve(cfdh[3][cfdh$plot==i,]*x+cfdh[2][cfdh$plot==i,],col=i,xlim=c(min(pine[pine$plot==i,]$ln_d2h),max(pine[pine$plot==i,]$ln_d2h)),add=TRUE)
	}

par(mfrow=c(2,4))
plot((lm(ln_Mst~ln_dbh+slope,data=pine)))
plot(lm(ln_Mst~ln_d2h+plot,data=pine))

summary(lm(ln_h~ln_dbh+slope,data=pine))


plot(summary(lm(ln_Mst~ln_d2h,data=pine))$res~pine$slope)
lines(lowess(pine$slope,summary(lm(ln_Mst~ln_d2h,data=pine))$res),col="blue")
plot(summary(lm(ln_h~ln_dbh,data=pine))$res~pine$std)


# quartz(w=12,h=8)
# par(mfrow=c(3,6))
# for (i in 1:max(as.numeric(bm$sp))){
	# plot(ln_h~ln_dbh,data=bm[bm$n.sp==i,],ylim=c(1,3.3),xlim=c(-4.1,-0.5))
	# title(bm$sp[bm$n.sp==i][1])
	# abline(lm(ln_h~ln_dbh,data=bm[bm$n.sp==i,]))
	# }


plot(ln_Mst~ln_d2h,data=pine)
lmo<-lm(ln_Mst~ln_d2h,data=pine)

plot(lmo)
abline(lmo)



