rm(list=ls())
getwd()
setwd('/Users/hyli0001/wrd/b/Dynamic_allometrics/')
# system('ls -alt ../processed_data')
bp<-readRDS('processed_data/plot_biomass.rds')
nfi_h<-readRDS('processed_data/nfi_harvested_sp.rds')
bp$sdi.ha<-2.4711*bp$sdi
bp$sdi_max.ha<-2.4711*bp$sdi_max
nfi_h$sdi.ha<-2.4711*nfi_h$sdi
bp[bp$sdi.ha>2000,]
## The number of spcies
sum(!is.na(unique(bp$sp_code)))

hist(nfi_h[nfi_h$sp_code=='CJ'&is.finite(nfi_h$sdi),'sdi.ha'],breaks=100)
nfi_h[nfi_h$sp_code=='CJ'&is.finite(nfi_h$sdi)&nfi_h$sdi.ha<10,]


## sdi.ha distribution

bp$dat_typ<-'hvst'
nfi_h$dat_typ<-'nfi'

quartz(w=8.0,h=5.2)
par(mfrow=c(2,1))
par(lwd=.3)
par(mai=c(.4,.5,.4,.1))
plot(NA, xlim=c(0.5,8.5),ylim=c(0,2500),xaxt='n',yaxt='n',xlab='',ylab='')
for (i in c('cnf')){
	n <- 1
	fml<-unique(bp[bp$ftp==i,'Family'])
	fml<-fml[order(fml)]
	for (j in fml){
		gns<-unique(bp[bp$ftp==i&bp$Family==j,'Genus'])
		gns<-gns[order(gns)]
		for (k in gns){
			spc<-unique(bp[bp$ftp==i&bp$Family==j&bp$Genus==k,'sp_code'])
			spc<-spc[order(spc)]
			for (l in spc){
				df_h<-bp[bp$sp_code==l,]
				df_nfi<-nfi_h[nfi_h$d.qm>0.001&nfi_h$sp_code==l,]
				df<-rbind(df_h[,c('ftp','sp_code','sdi.ha','dat_typ')]
,df_nfi[,c('ftp','sp_code','sdi.ha','dat_typ')])

df<-df[is.finite(df$sdi.ha),]
df$x<-runif(length(df$sdi.ha),min=-0.1,max=0.1)+n
y.m<-mean(df$sdi.ha)
y.sd<-sd(df$sdi.ha)
y.lim<-range(df$sdi.ha)
y.n<-length(df$sdi.ha)

y.grid<-seq(y.lim[1],y.lim[2],length.out=250)
y.prop.dist<-dnorm(y.grid,y.m,y.sd)
y.norm<-y.prop.dist/(max(y.prop.dist)*5)
sdi_max<-unique(df_h$sdi_max.ha)

# data points
points(sdi.ha~x,df[df$dat_typ=='nfi',],xlim=c(-1,3),lwd=0.25,cex=0.3,col=8)
points(sdi.ha~x,df[df$dat_typ=='hvst',],xlim=c(-1,3),lwd=0.25,cex=0.3,pch=21,bg=8)

# boxplots
boxplot(df$sdi.ha,at=n,add=TRUE,boxwex=0.75,outline=FALSE,axes=FALSE,col=0,border="black",lwd=1)

### 95th percentile -- sdi_max
points(n,sdi_max,lwd=1.5,lty=2,pch=21,bg=2)

# ### violin
# if(length(df$sdi.ha) > 1 && length(unique(df$sdi.ha)) > 1){
	# y.den <- density(df$sdi.ha)
	# y.grid <- y.den$x
	# y.norm <- y.den$y / max(y.den$y) * 0.20
# lines(y.norm+n, y.grid)
# lines(-y.norm+n, y.grid)
# }

n <- n+1
}}}}

par(lwd=.3)
par(mai=c(.4,.5,.4,.1))
plot(NA, xlim=c(0.5,14.5),ylim=c(0,2500),xaxt='n',yaxt='n',xlab='',ylab='')
for (i in c('brd')){
	n <- 1
	fml<-unique(bp[bp$ftp==i,'Family'])
	fml<-fml[order(fml)]
	for (j in fml){
		gns<-unique(bp[bp$ftp==i&bp$Family==j,'Genus'])
		gns<-gns[order(gns)]
		for (k in gns){
			spc<-unique(bp[bp$ftp==i&bp$Family==j&bp$Genus==k,'sp_code'])
			spc<-spc[order(spc)]
			for (l in spc){
				df_h<-bp[bp$sp_code==l,]
				df_nfi<-nfi_h[nfi_h$d.qm>0.001&nfi_h$sp_code==l,]
				df<-rbind(df_h[,c('ftp','sp_code','sdi.ha','dat_typ')]
,df_nfi[,c('ftp','sp_code','sdi.ha','dat_typ')])

df<-df[is.finite(df$sdi.ha),]
df$x<-runif(length(df$sdi.ha),min=-0.1,max=0.1)+n
y.m<-mean(df$sdi.ha)
y.sd<-sd(df$sdi.ha)
y.lim<-range(df$sdi.ha)
y.n<-length(df$sdi.ha)

y.grid<-seq(y.lim[1],y.lim[2],length.out=250)
y.prop.dist<-dnorm(y.grid,y.m,y.sd)
y.norm<-y.prop.dist/(max(y.prop.dist)*5)
sdi_max<-unique(df_h$sdi_max.ha)

# data points
points(sdi.ha~x,df[df$dat_typ=='nfi',],xlim=c(-1,3),lwd=0.25,cex=0.3,col=8)
points(sdi.ha~x,df[df$dat_typ=='hvst',],xlim=c(-1,3),lwd=0.25,cex=0.3,pch=21,bg=8)

# boxplots
boxplot(df$sdi.ha,at=n,add=TRUE,boxwex=0.75,outline=FALSE,axes=FALSE,col=0,border="black",lwd=1)

### 95th percentile -- sdi_max
points(n,sdi_max,lwd=1.5,lty=2,pch=21,bg=2)

# ### violin
# if(length(df$sdi.ha) > 1 && length(unique(df$sdi.ha)) > 1){
	# y.den <- density(df$sdi.ha)
	# y.grid <- y.den$x
	# y.norm <- y.den$y / max(y.den$y) * 0.20
# lines(y.norm+n, y.grid)
# lines(-y.norm+n, y.grid)
# }

n <- n+1
}}}}



	
unique(bp$Family)
unique(bp$Genus)
unique(bp$sp_code)

for (i in 1:22){
	df<-nfi_h[nfi_h$d.qm>0.001&nfi_h$sp_code==unique(nfi_h$sp_code)[i],]
	
	
	

y<-df$sdi[is.finite(df$sdi)]
y.m<-mean(y)
y.sd<-sd(y)
y.lim<-range(y)
y.n<-length(y)

x<-runif(y.n,min=-0.1,max=0.1)
y.grid<-seq(y.lim[1],y.lim[2],length.out=250)
y.prop.dist<-dnorm(y.grid,y.m,y.sd)
y.norm<-y.prop.dist/(max(y.prop.dist)*5)

df_h<-bp[bp$sp_code==unique(nfi_h$sp_code)[12],]
y_h<-df_h$sdi[is.finite(df_h$sdi)]
x_h<-runif(length(y_h),min=-0.1,max=0.1)

plot(x,y,xlim=c(-1,3),col=8,lwd=0.1,cex=0.3)
points(x_h,y_h,xlim=c(-1,3),lwd=0.25,cex=0.3,pch=21,bg=8)
lines(y.norm,y.grid)
lines(-y.norm,y.grid)


head(dnorm(y.grid,y.m,y.sd))
norm.density<-dnorm(y_grid, sdi.m, sdi.sd)/max(dnorm(y_grid, sdi.m, sdi.sd))

plot(y_grid, norm.density,ylim=c(0,1),xlim=c(0,800))

## Gamma distribution : shape=mean^2/var; rate=mean/var; scale=var/mean; Dispersal=var/mean^2

quartz(w=4.5,h=2.55)
par(mfrow=c(1,4))
par(lwd=.1,col=0)
par(mai=c(.35,0,.2,.0))
plot(NA,xlab="",ylab="",yaxt="n",xaxt="n",xlim=c(0,10),ylim=c(-.2,9.15))
# ,ylim=c(.5,8.35)
x.label<-unique(smg[,c('spe','rev.ord')])
x.label<-x.label[order(-x.label$rev.ord),]
x.label$ful.nm<-c('Pinus sylvestris','Pinus taeda','Eucalyptus grandis','Pseudotsuga menziesii','Populus tremula × tremuloides','Betula pendula','Bruguiera gymnorrhiza','Picea abies','Broad-leaved deciduous')

par(lwd=.1,col=1)
mtext(x.label$ful.nm[1:8],las=2,2,at=(x.label$rev.ord[1:8])-0.1,line=-8.5,font=3,cex=6/12)
mtext(x.label$ful.nm[9],las=2,2,at=(x.label$rev.ord[9])-0.1,line=-8.5,font=c(1),cex=6/12)
mtext('Average',las=2,2,at=-0.25,line=-8,font=c(1),cex=6/12)

for (i in c(6:9)){
	lnb<-rbind(lnd[1,],lnd[lnd$rev.ord%in%i,])
	lnb[1,]<-NA
legend(yjust=1.25,3,i,x.intersp=.4,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}
for (i in c(3:5)){
	lnb<-rbind(lnd[1:3,],lnd[lnd$rev.ord%in%i,])
lnb[1:3,]<-NA
legend(yjust=1.25,3,i,x.intersp=.4,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}

for (i in c(1)){
	lnb<-rbind(lnd[1:3,],lnd[lnd$rev.ord%in%i,])
lnb[1:3,]<-NA
legend(yjust=1.25,3,i,x.intersp=.4,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}

for (i in c(2)){
lnb<-lnd[lnd$rev.ord%in%i,][1:7,]
legend(yjust=1.2,-.3,i,x.intersp=.4,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=.8)
}


par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,0.9),ylim=c(-.2,9.25))
points((jitter(rev.ord,.5)-0.5)~brtr,smg,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(smg$sp==1,bg,col))
mtext('a',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-4,4,by=0.5),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('∆'[B],' (y'^-1,')')),1,line=.9,font=1,cex=7/12)

for (i in 1:2){
	df<-smg[smg$sp==i,]
	for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$brtr[dff$brtr>=0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,2,0.001),shape,rate))*1.85
	y.mn<-unique(dff$rev.ord)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
	# print(unique(dff[,c('site','sp')]))
	# print(gamma_test(x))
}}

	for (i in 1:nrow(smm)){
		lines(c(-1,1)*smm$brtr.sd[i]+smm$brtr.m[i],rep(smm$rev.ord[i],2)-.5)}
	points((rev.ord-.5)~brtr.m,smm,bg=bg,pch=21,col=col,cex=7/12,lwd=ifelse(smm$sp==1,.2,.7))













year<-1991:2002
Tb.dbc<-data.frame(year)
Tb.dbc$a<-NA
Tb.dbc$b<-NA
Tb.dbc$p<-NA
Tb.dbc$lc.m<-NA
Tb.dbc$lc.se<-NA
Tb.dbc$bc.d.m<-NA
Tb.dbc$bc.d.se<-NA
for (i in 1991:2002){
	df<-md.ob.in[md.ob.in$year==i&md.ob.in$bc.d>0.01,]
	coe0<-summary(lm(bc.d~lc+0,df))$coe
	Tb.dbc[Tb.dbc$year==i,c('b0')]<-coe0[1,1]
	Tb.dbc[Tb.dbc$year==i,c('p0')]<-coe0[1,4]

	coe<-summary(lm(bc.d~lc,df))$coe
	Tb.dbc[Tb.dbc$year==i,c('a','b')]<-coe[1:2,1]
	Tb.dbc[Tb.dbc$year==i,c('a.p')]<-coe[1,4]
	Tb.dbc[Tb.dbc$year==i,c('b.p')]<-coe[2,4]
	Tb.dbc[Tb.dbc$year==i,c('lc.m','lc.se','bc.d.m','bc.d.se')]<-summaryBy(lc+bc.d~year,df,FUN=me)[,c('lc.m','lc.se','bc.d.m','bc.d.se')]
	Tb.dbc[Tb.dbc$year==i,c('lc.mn','lc.mx')]<-summaryBy(lc~year,df,FUN=range)[2:3]
	}
Tb.dbc$func.0<-ifelse(Tb.dbc$p0<0.05,1,0)	
Tb.dbc$func.1<-ifelse(Tb.dbc$a.p<0.05,1,0)
Tb.dbc$func.2<-ifelse(Tb.dbc$b.p<0.05,1,0)


# Modelling evaluation using Tim's data (SETRES)
anova(lme(p.Bbr~as.factor(year):as.factor(strata)+as.factor(plot)+ as.factor(year)*as.factor(strata),dis.inv.str,random=~1|block))

anova(lme(p.Bdb~as.factor(year):as.factor(strata)+as.factor(plot)+ as.factor(year)*as.factor(strata),dis.inv.str,random=~1|block))

quartz(w=3.42,h=1.8)
par(mfrow=c(1,2))
par(lwd=.3)
year<-1991:2002

par(mai=c(.4,.4,.2,.02))
plot(NA,xlim=c(0,.27),ylim=c(0,1.1),xlab='',ylab='',xaxt='n',yaxt='n')
mtext('a',line=-.8,adj=0.02,font=2,cex=8/12)
mtext('Pinus taeda',line=-.1,adj=0,font=3,cex=7/12)

for (i in year) for (j in 1:4){
	df<-dis.Bdb.y[dis.Bdb.y$year==i&dis.Bdb.y$plot==j,]
	points(rlcn~p.Bbr.m,df,col=j,pch=21,cex=.4,lwd=.05,type='o')}

df<-dis.Bdb	
	points(rlcn~p.Bbr.m.mean,summaryBy(p.Bbr.m~rlcn,df,FUN=mean),type='o',lwd=.5,cex=0)
	for (j in 1:nrow(df)){
		lines(df$p.Bbr.m[j]+c(-1,1)*df$p.Bbr.se[j],rep(df$rlcn[j],2))}
		points(rlcn~p.Bbr.m,df,bg=plot,pch=21,cex=7/12,lwd=.2)
axis (1,seq(-1,2,by=.1),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=.3)
axis (2,seq(-1,2,by=.5),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=.3)
mtext(expression(paste('Proportional branch biomass')),side=1,line=.5,font=1,cex=7/12)
mtext(expression(paste('Relative branch position')),side=2,line=1,font=1,cex=7/12)
mtext(expression(paste('within crown')),side=2,line=.6,font=1,cex=7/12)

par(mai=c(.4,.02,.2,.4))
plot(NA,xlim=c(0,1),ylim=c(0,1.1),col=0,xlab='',ylab='',xaxt='n',yaxt='n')
mtext('b',line=-.8,adj=0.02,font=2,cex=8/12)
mtext('SETRES, 1991-2002',line=-.1,adj=1,font=1,cex=7/12)

for (i in year) for (j in 1:4){
	df<-dis.Bdb.y[dis.Bdb.y$year==i&dis.Bdb.y$plot==j,]
	points(rlcn~p.Bdb.m,df,col=j,pch=21,cex=.4,lwd=.05,type='o')}

df<-dis.Bdb
	points(rlcn~p.Bdb.m.mean,summaryBy(p.Bdb.m~rlcn,df,FUN=mean),type='o',lwd=.5,cex=0)
	for (j in 1:nrow(df)){
		lines(df$p.Bdb.m[j]+c(-1,1)*df$p.Bdb.se[j],rep(df$rlcn[j],2))}
		points(rlcn~p.Bdb.m, df,bg=plot,pch=21,cex=7/12,lwd=.2)
axis (1,seq(-1,2,by=.5),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=.3)
axis (2,seq(-1,2,by=.5),tck=.02,label=FALSE,mgp=c(0,0,0),cex.axis=7/12,lwd=.3)
mtext(expression(paste('Proportional branch turnover')),side=1,line=.5,font=1,cex=7/12)

legend('topright',x.intersp=.5,y.intersp=.85,c('C','I','F','IF'),pch=21,pt.bg=c(1:4),col=1,box.col=0,cex=7/12,pt.lwd=.3,horiz=TRUE)

n.sta<-seq(3,8,by=1)
n.col<-8
itv<-1
dff<-mo.id
dff[,c('brt.kg','Bdb.kg.sum')]<-mo.id[,c('brt.g','Bdb.g.sum')]/1000
dff<-dff[order(dff$brt.g),]
df<-dff[dff$var=='tBbr.g'&dff$n.sta%in%n.sta&dff$itv==itv,]
crv<-stv[stv$var=='tBbr.g'&stv$str%in%n.sta&stv$itv==itv,]
quartz(w=3.42,h=1.8)
par(mfrow=c(1,2))
par(lwd=.3)

par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,12.5),xlim=c(0,12.5))
abline(0,b=1,lty=3,lwd=.75)
mtext('a',line=-.8,adj=0.02,font=2,cex=8/12)
mtext(expression(paste('Gross biomass')),line=-.8,adj=0.5,font=1,cex=7/12)
mtext('Pinus taeda',line=-.1,adj=0,font=3,cex=7/12)

points(brt.kg~Bdb.kg.sum,df,pch=1,col=terrain.colors(n.col)[n.sta],cex=4/12,lwd=.05)
for(i in n.sta){
	curve(crv[crv$str==i,'slp']*x,add=TRUE,xlim=range(df$Bdb.kg.sum,na.rm=TRUE),col=terrain.colors(n.col)[i],lwd=1)
	}
points(brt.kg~Bdb.kg.sum,summaryBy(brt.kg+Bdb.kg.sum~plot+block+Tnr,df,FUN=mean,keep.names=TRUE),pch=21,bg=terrain.colors(n.col)[n.sta],cex=7/12,lwd=.2)
axis (1,seq(-3,15,by=3),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=.3)
axis (2,seq(-3,15,by=3),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=.3)
mtext(expression(paste('Turnover'['B'][' mod'],' (kg tree'^-1,' y'^-1,')')),side=2,line=.5,font=1,cex=7/12)


par(mai=c(.575,1.23,.85,.04),new=TRUE)
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(2,9),ylim=c(0.5,1.2))
lines(c(0,20),c(1,1),lty=1)
points(slp~str,stv[stv$str%in%n.sta&stv$var=='tBbr.g',],cex=4.5/12,type='o',pch=21,bg=terrain.colors(n.col)[str],lwd=.2)
axis (1,seq(-10,20,by=3),tck=.02,label=TRUE,mgp=c(0,-.4,0),cex.axis=6/12,lwd=.3)
axis (2,seq(.5,1.2,by=.5),tck=.02,label=TRUE,mgp=c(0,-.1,0),cex.axis=6/12,lwd=.3)
mtext(expression(paste('strata')),side=1,line=-.1,font=1,cex=6/12)
mtext(expression(paste('slope')),side=2,line=.3,font=1,cex=6/12)

df<-dff[dff$var=='Bbr.g'&dff$n.sta%in%n.sta&dff$itv==itv,]
crv<-stv[stv$var=='Bbr.g'&stv$str%in%n.sta&stv$itv==itv,]

par(mai=c(.4,.02,.2,.4))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,12.5),xlim=c(0,12.5))
abline(0,b=1,lty=3,lwd=.75)
mtext('b',line=-.8,adj=0.02,font=2,cex=8/12)
mtext(expression(paste('Live biomass')),line=-.8,adj=0.5,font=1,cex=7/12)
mtext('SETRES, 1991-2002',line=-.1,adj=1,font=1,cex=6/12)

points(brt.kg~Bdb.kg.sum,df,pch=1,col=terrain.colors(n.col)[n.sta],cex=4/12,lwd=.05)
for(i in n.sta){
	curve(crv[crv$str==i,'slp']*x,add=TRUE,xlim=range(df$Bdb.kg.sum,na.rm=TRUE),col=terrain.colors(n.col)[i],lwd=.5)
	}
points(brt.kg~Bdb.kg.sum,summaryBy(brt.kg+Bdb.kg.sum~plot+block+Tnr,df,FUN=mean,keep.names=TRUE),pch=21,bg=terrain.colors(n.col)[n.sta],cex=7/12,lwd=.2)
axis (1,seq(-3,15,by=3),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=.3)
axis (2,seq(-3,15,by=3),tck=.02,label=FALSE,mgp=c(0,0,0),cex.axis=7/12,lwd=.3)

mtext(expression(paste('Observed annual branch turnover (kg tree'^-1,' y'^-1,')')),side=1,line=.6,font=1,at=-1,cex=7/12)

par(mai=c(.575,.85,.85,.42),new=TRUE)
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(2,9),ylim=c(0.5,1.2))
lines(c(0,20),c(1,1),lty=1)
points(slp~str,stv[stv$str%in%n.sta&stv$var=='Bbr.g',],cex=4.5/12,type='o',pch=21,bg=terrain.colors(n.col)[str],lwd=.2)
axis (1,seq(-10,20,by=3),tck=.02,label=TRUE,mgp=c(0,-.4,0),cex.axis=6/12,lwd=.3)
axis (2,seq(.5,1.2,by=.5),tck=.02,label=TRUE,mgp=c(0,-.1,0),cex.axis=6/12,lwd=.3)
mtext(expression(paste('strata')),side=1,line=-.1,font=1,cex=6/12)
mtext(expression(paste('slope')),side=2,line=.3,font=1,cex=6/12,at=0.8)



sens<-mo.id[mo.id$var=='Bbr.g'&mo.id$n.sta==3&mo.id$itv==1,]
sens[,c('brt.kg','Bdb.kg.sum')]<-mo.id[mo.id$var=='Bbr.g'&mo.id$n.sta==3&mo.id$itv==1,c('brt.g','Bdb.g.sum')]/1000
sens$Bbr.p<-sens$Bbr.g.m/sens$Bbr.g.sum
anova(lme(brt.kg~Bdb.kg.sum+Bbr.g.sum,sens,random=~1|block))

# No effect on the % within-canopy mortality of treatment/year/size


## Population density of Biomass density Using measured individuals
tnr<-Tnr.nd[Tnr.nd$site!='TECHS29',]
tnr$rev.ord<-(10)-tnr$ord.spe
tnr$rev.ord.1<-tnr$rev.ord
tnr$rev.ord.1[tnr$rev.ord>3]<-tnr$rev.ord[tnr$rev.ord>3]-1

tnr$rc<-tnr$lc/tnr$h
tnm<-summaryBy(Bbr.kg.m+bc.d+rc+brtr+h.d~bg+col+ord.spe+rev.ord.1+sp+spe,tnr,FUN=md)
tnm$ord.spe<-jitter(tnm$ord.spe,ifelse(tnm$ord.spe==6,0,.5))
tnm$rev.ord.1<-jitter(tnm$rev.ord.1,ifelse(tnm$rev.ord.1==1,0,.5))

lnd<-unique(lnd[,c('sp','spe','site','bg','col','ord','ord.spe','bg.r','col.r','lwd','lty','st')])
lnd$rev.ord<-10-lnd$ord.spe
lnd$rev.ord.1<-lnd$rev.ord
lnd$rev.ord.1[lnd$rev.ord>3]<-lnd$rev.ord[lnd$rev.ord>3]-1

# write.table(tnr,'/Users/hyli0001/Documents/wd/7_Branch_turnover/BT_Ind_values_for_non.linear.dat',sep='\t',quote=FALSE,row.names=FALSE)

## Gamma distribution : shape=mean^2/var; rate=mean/var; scale=var/mean; Dispersal=var/mean^2
library('fitdistrplus')
library('goft')
# par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# fw <- fitdist(dff$Bbr.kg.m[dff$Bbr.kg.m>0], "weibull")
# fln <- fitdist(dff$Bbr.kg.m[dff$Bbr.kg.m>0], "lnorm")
# fg <- fitdist(dff$Bbr.kg.m[dff$Bbr.kg.m>0], "gamma")
# plot.legend <- c("Weibull", "lognormal", "gamma")
# denscomp(list(fw, fln, fg), legendtext = plot.legend)
# qqcomp(list(fw, fln, fg), legendtext = plot.legend)
# cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
# ppcomp(list(fw, fln, fg), legendtext = plot.legend)


quartz(w=4.5,h=2.2)
par(mfrow=c(1,4))
par(lwd=.1,col=0)
par(mai=c(.35,0,.2,.0))
plot(NA,xlab="",ylab="",yaxt="n",xaxt="n",xlim=c(0,10),ylim=c(.5,8.35))

x.label<-unique(tnr[,c('spe','rev.ord')])
x.label$rev.ord.1<-x.label$rev.ord
x.label$rev.ord[x.label$rev.ord.1>3]<-x.label$rev.ord.1[x.label$rev.ord.1>3]-1
x.label<-x.label[order(-x.label$rev.ord),]
x.label$ful.nm<-c('Pinus sylvestris','Pinus taeda','Eucalyptus grandis','Pseudotsuga menziesii','Populus tremula × tremuloides','Bruguiera gymnorrhiza','Picea abies','Broad-leaved deciduous')
par(lwd=.1,col=1)
mtext(x.label$ful.nm[1:7],las=2,2,at=(x.label$rev.ord[1:7])-0.1,line=-8.5,font=3,cex=6/12)
mtext(x.label$ful.nm[8],las=2,2,at=(x.label$rev.ord[8])-0.1,line=-8.5,font=c(1),cex=6/12)
lnd$rev.ord.1<-lnd$rev.ord
lnd$rev.ord.1[lnd$rev.ord>3]<-lnd$rev.ord[lnd$rev.ord>3]-1
lnd$rev.ord.1[lnd$site=='ja']<-NA
rbind(lnd[1,],lnd[lnd$rev.ord.1%in%7,])
for (i in c(5:8)){
	lnb<-rbind(lnd[1,],lnd[lnd$rev.ord.1%in%i,])
	lnb[1,]<-NA
legend(yjust=1.25,3.5,i,x.intersp=.3,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}
for (i in c(3:4)){
	lnb<-rbind(lnd[1:3,],lnd[lnd$rev.ord.1%in%i,])
lnb[1:3,]<-NA
legend(yjust=1.25,3.5,i,x.intersp=.3,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}

for (i in c(1)){
	lnb<-rbind(lnd[1:3,],lnd[lnd$rev.ord.1%in%i,])
lnb[1:3,]<-NA
legend(yjust=1.25,3.5,i,x.intersp=.3,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}

for (i in c(2)){
lnb<-lnd[lnd$rev.ord.1%in%i,][1:7,]
legend(yjust=1.25,-.3,i,x.intersp=.3,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=.8)
}

par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,3.5),ylim=c(.5,8.35))
points(jitter(rev.ord.1,.5)-.5~brtr,tnr,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(tnr$sp==1,bg,col))
mtext('a',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-4,15,by=1),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('†'[BR],' (yr'^-1,')')),1,line=.9,font=1,cex=7/12)

for (i in 1:2){
df<-tnr[tnr$sp==i,]
for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$brtr[dff$brtr>0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-ifelse(dff$site=='ok',c(0,2),range(x,na.rm=TRUE))
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,2,0.001),shape,rate))*1.75
	y.mn<-unique(dff$rev.ord.1)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.3),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
	}}
for (i in 1:nrow(tnm)){
		lines(c(-1,1)*tnm$brtr.sd[i]+tnm$brtr.m[i],rep(tnm$rev.ord.1[i],2)-.5)}
	points((rev.ord.1-.5)~brtr.m,tnm,bg=bg,pch=21,col=col,cex=7/12,lwd=ifelse(tnm$sp==1,.2,.7))

par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,11),ylim=c(.5,8.35))
points(jitter(rev.ord.1,.5)-.5~bc.d,tnr,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(tnr$sp==1,bg,col))
mtext('b',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-3,15,by=3),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('Crown ascent rate')),1,line=.7,font=1,cex=7/12)
mtext(expression(paste(' (m y'^-1,')')),1,line=1.6,font=1,cex=7/12)

for (i in 1:2){
df<-tnr[tnr$sp==i,]
for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$bc.d[dff$bc.d>0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,10,0.001),shape,rate))*1.75
	y.mn<-unique(dff$rev.ord.1)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
	}}

	for (i in 1:nrow(tnm)){
		lines(c(-1,1)*tnm$bc.d.sd[i]+tnm$bc.d.m[i],rep(tnm$rev.ord.1[i],2)-.5)}
	points((rev.ord.1-.5)~bc.d.m,tnm,bg=bg,pch=21,col=col,cex=7/12,lwd=ifelse(tnm$sp==1,.2,.7))

par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,8),ylim=c(.5,8.35))
points((jitter(rev.ord.1,.5)-0.5)~Bbr.kg.m,tnr,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(tnr$sp==1,bg,col))
mtext('c',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-4,15,by=2),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('Branch biomass density')),1,line=.7,font=1,cex=7/12)
mtext(expression(paste('(kg m'^-1,')')),1,line=1.6,font=1,cex=7/12)
mtext('Measured individual trees',line=-.1,adj=1,font=1,cex=7/12)

for (i in 1:2){
df<-tnr[tnr$sp==i,]
for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$Bbr.kg.m[dff$Bbr.kg.m>0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,10,0.001),shape,rate))*1.75
	y.mn<-unique(dff$rev.ord.1)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)

# # 	print(unique(dff[,c('site','sp')]))
	# print(gamma_test((x)))
	}}

	for (i in 1:nrow(tnm)){
		lines(c(-1,1)*tnm$Bbr.kg.m.sd[i]+tnm$Bbr.kg.m.m[i],rep(tnm$rev.ord.1[i],2)-.5)}
	points((rev.ord.1-.5)~Bbr.kg.m.m,tnm,bg=bg,pch=21,col=col,cex=7/12,lwd=ifelse(tnm$sp==1,.2,.7))


lnd0<-lnd[1,]
lnd0[1,]<-NA
lnd5<-lnd[1:5,]
lnd5[1:5,]<-NA
lnd4<-lnd[1:4,]
lnd4[1:4,]<-NA
lnd3<-lnd[1:3,]
lnd3[1:3,]<-NA
lnd6<-lnd[1:6,]
lnd6[1:6,]<-NA


# # Analyzing factors explaining variation of the turnover 

btdk1<-btd[btd$site=='dk'&btd$sp==1,c('trt1','trt2','plot','std')]
btdk2r<-btd[btd$site=='dk'&btd$sp==2,c('trt1','trt2','plot','std','site','sp')]
btdk2<-merge(btdk1,btdk2r,by=c('trt1','trt2','plot'))
btdk2$std2<-btdk2$std.x+btdk2$std.y
btdk2$std.x<-NULL
btdk2$std.y<-NULL

btfr1<-btd[btd$site=='fr'&btd$trt1=='M'&btd$sp==1,c('trt1','plot','std')]
btfr2r<-btd[btd$site=='fr'&btd$trt1=='M'&btd$sp==2,c('trt1','plot','std','site','trt2','sp')]
btfr2<-merge(btfr1,btfr2r,by=c('trt1','plot'))
btfr2$std2<-btfr2$std.x+btfr2$std.y
btfr2$std.x<-NULL
btfr2$std.y<-NULL

smg.spr<-rbind(btfr2, btdk2)
smg<-merge(btd,smg.spr,by=c('trt1','trt2','plot','site','sp'),all.x=TRUE)
smg$std.sp<-smg$std
smg$std<-ifelse(!is.na(smg$std2),smg$std2,smg$std.sp)
smg$ttrt<-paste(smg$trt1,smg$trt2,sep='')
smg$std.mx<-exp(log((smg$Bbr*10/smg$std)/(2.528*1.1))*(1/-1.515))*10000
smg$sdi<-(smg$std/2.4711)*(smg$d*100/25.4)^ifelse(smg$spe%in%c('P. abies','P. sylvestris','P. tremula × tremuloides','B. pendula'),1.66,ifelse(smg$spe=='B. gymnorrhiza',1.6175,1.605))

smg$rsd1<-ifelse(smg$spe=='B. gymnorrhiza',smg$std,smg$sdi)/ifelse(smg$spe=='B. gymnorrhiza',smg$std.mx,ifelse(smg$spe=='P. taeda',450,400))
smg$rsd<-ifelse(smg$spe=='B. gymnorrhiza',smg$std,smg$sdi)/ifelse(smg$spe=='B. gymnorrhiza',smg$std.mx,400)

smg$NPPbr<-smg$NPPbr0+smg$brt
smg$NPPaw<-smg$NPPaw0+smg$brt
smg$brtr1<-smg$brt/smg$Bbr.y
smg$brtr.nl<-smg$brt/smg$Bbr.y
smg$h.d.lc<-smg$h.d/smg$lc
smg$rc<-smg$lc/smg$h
smg$bc.d.lc<-smg$bc.d/smg$lc
smg$dh<-smg$d/smg$h
smg$rh<-smg$h.d/smg$h
smg$rNPPbr<-smg$brt/smg$NPPbr
smg$rNPPaw<-smg$brt/smg$NPPaw
smg$rNPPaw0<-smg$NPPaw0/smg$NPPaw
smg$rNPPaw01<-smg$NPPaw/smg$NPPaw0
smg$rev.ord<-(10)-smg$ord.spe
smm<-summaryBy(rNPPbr+rNPPaw+rc+brtr+h.d~bg+col+ord.spe+rev.ord+sp+spe,smg,FUN=md)
smm$ord.spe<-jitter(smm$ord.spe,ifelse(smm$ord.spe==6,0,.5))
smm$rev.ord<-jitter(smm$rev.ord,ifelse(smm$rev.ord==1,0,.5))

smg$Bbr.0<-smg$Bbr.y-smg$brt
smg$brt.within.crown<-(smg$brt*0.16/0.84)

mean(smg$brt.within.crown/smg$NPPaw0,na.rm=TRUE)
sd(smg$brt.within.crown/smg$NPPaw0,na.rm=TRUE)

sd(smg$brtr)/mean(smg$brtr)
mean(smg$brtr[smg$site=='sl'])
sd(smg$brtr[smg$site=='sl'])
range(smg$rNPPaw[smg$site=='sl'])
md(smg$brtr)
range(smg$brtr)
range(smg$rNPPbr)
range(smg$rNPPaw)
md(smg$rNPPbr)
md(smg$rNPPaw)
md(smg$rNPPaw01)
md(smg$rNPPaw01[smg$site%in%c('TECHS22','TECHS20','TECHS30','ok')])
md(smg$rNPPaw01[smg$site%in%c('dk','fv','ml','mt','sl','nc')])
md(smg$rNPPaw01[smg$site%in%c('hb','ro','Bräc','Gävl','Grän','Möln','Ebbe','ja','fl','fr','rd')])

range(smg$brtr[smg$sp==1])
range(smg$brtr[smg$sp==2])
mean((smg[smg$site=='dk'&smg$sp==1,][1:4,'brt']/0.84)/smg[smg$site=='dk'&smg$sp==1,][1:4,'NPPaw'])

write.table(smg,'/Users/hyli0001/Documents/wd/7_Branch_turnover/BT_Stand_values_for_non.linear.dat',sep='\t',quote=FALSE,row.names=FALSE)


######################################################
## Logistic prediction ##############################
h.d<-function(h.d,sp){
	mapply(function(h.d,sp){
		if(sp==1){
			brtr.h.d<-nls(brtr~c+((d-c)/(1+exp(-a*(h.d-b)))),smg[smg$sp==1,],start=list(a=1.74,b=2.55,c=.058,d=.73))
				a<-summary(brtr.h.d)$coe[1]
				b<-summary(brtr.h.d)$coe[2]
				c<-summary(brtr.h.d)$coe[3]
				d<-summary(brtr.h.d)$coe[4]
			return(c+((d-c)/(1+exp(-a*(h.d-b)))))
				}
		if(sp==2){
		brtr.h.d<-lm(brtr~h.d,smg[smg$sp==2,])
				coe<-summary(brtr.h.d)$coe
			return(coe[1]+coe[2]*h.d)
				}
							},h.d,sp)
					}


## Prediction and residuals
smg$h.d.pre<-h.d(smg$h.d,smg$sp)
smg$h.d.res<-smg$brtr-smg$h.d.pre
smg$h.d.res.r<-smg$brtr/smg$h.d.pre
summary(lm(h.d.pre~brtr,smg[smg$sp==1,]))
## Relative residulas vs. std
h.d.std.r<-function(std,rsd,sp,nl){
	mapply(function(std,rsd,sp,nl){
	if(sp==1){
		if(nl==1){
		coe<-
		nls(h.d.res.r~c+((d-c)/(1+exp(-a*(rsd-b)))),smg[smg$sp==1,],start=list(a=6.06,b=.7212,c=.5796,d=1.6384))
			a<-summary(coe)$coe[1]
			b<-summary(coe)$coe[2]
			c<-summary(coe)$coe[3]
			d<-summary(coe)$coe[4]
			return(c+((d-c)/(1+exp(-a*(rsd-b)))))
				}
		if(nl==2){
		coe<-summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))$coe
		return(coe[1]+(rsd)*coe[2])
				}
				}
	if(sp==2){
	coe<-summary(lm(h.d.res.r~std,smg[smg$sp==2,]))$coe
	return(coe[1]+(std)*coe[2])}	
	},std,rsd,sp,nl)
	}

smg$rsd.pre<-h.d.std.r(smg$std,smg$rsd,smg$sp,2)
summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))
smg$rsd.res<-smg$h.d.res.r-smg$rsd.pre


# Residual distribution

########## Whole model ################
h.d.std<-function(h.d,sp,std,rsd,nl){
	mapply(function(h.d,sp,std,rsd,nl){
		if(sp==1){
			brtr.h.d<-nls(brtr~c+((d-c)/(1+exp(-a*(h.d-b)))),smg[smg$sp==1,],start=list(a=2.22,b=2.55,c=.09,d=.81))
				a<-summary(brtr.h.d)$coe[1]
				b<-summary(brtr.h.d)$coe[2]
				c<-summary(brtr.h.d)$coe[3]
				d<-summary(brtr.h.d)$coe[4]
		if(nl==1){
		coe<-
		nls(h.d.res.r~c+((d-c)/(1+exp(-a*(rsd-b)))),smg[smg$sp==1,],start=list(a=6.06,b=.7212,c=.5796,d=1.6384))
			a1<-summary(coe)$coe[1]
			b1<-summary(coe)$coe[2]
			c1<-summary(coe)$coe[3]
			d1<-summary(coe)$coe[4]
		return((c+((d-c)/(1+exp(-a*(h.d-b)))))*(c1+((d1-c1)/(1+exp(-a1*(rsd-b1))))))
				}
		if(nl==2){
		coe<-summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))$coe
		return((c+((d-c)/(1+exp(-a*(h.d-b)))))*(coe[1]+(rsd)*coe[2]))
				}
				}

		if(sp==2){
		brtr.h.d<-lm(brtr~h.d,smg[smg$sp==2,])
				coe1<-summary(brtr.h.d)$coe
		brtr.res.r.std2<-lm(h.d.res.r~std,smg[smg$sp==2,])
				coe2<-summary(brtr.res.r.std2)$coe
			return((coe1[1]+coe1[2]*h.d)*(coe2[1]+coe2[2]*std))
				}
							},h.d,sp,std,rsd,nl)
					}

smg$pre.brtr<-h.d.std(smg$h.d,smg$sp,smg$std,smg$rsd,2)
smg$res.brtr<-smg$brtr-smg$pre.brtr
smg$res.brtr.r<-smg$brtr/smg$pre.brtr
# # # # # TESTING! Treatment effects  
# # # # for (i in c("P. abies","E. grandis","P. taeda","P. sylvestris","P. menziesii")){
	# # df<-smg[smg$site==i&smg$res.brtr.r>0,]
	
	# # anv<-anova(glm(res.brtr.r~site,smg[smg$res.brtr.r>0,],family=Gamma(link='identity')))
	# # anv$p.v<-round((1-pf((anv[2,2]/anv[2,1])/(anv[2,4]/anv[2,3]), anv[2,1], anv[2,3])),4)}



quartz(w=3.42,h=1.8)
par(mfrow=c(1,2))
par(lwd=.3)
par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,1.1),ylim=c(-0.3,0.3))
lines(c(-1,2),c(0,0),col=8,lwd=7/12,lty=3)
points(res.brtr~pre.brtr,smg,col=col,pch=pch,bg=bg,cex=.5,lwd=.2)
points(res.brtr~pre.brtr,smg[smg$sp==2&smg$site=='dk',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)
points(res.brtr~pre.brtr,smg[smg$sp==2&smg$site=='fr',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)
points(res.brtr~pre.brtr,smg[smg$sp==1&smg$site=='ok',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)	
	
mtext('a',line=-.8,adj=0.02,font=2,cex=8/12)
axis (1,seq(-2,2,by=.5),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-2,2,by=.2),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Residuals')),2,line=.7,font=1,cex=7/12)
mtext(expression(paste('Predicted ‡'[B])),1,line=.5,font=1,cex=7/12)

par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,1.1),ylim=c(0,3.5))
lines(c(-1,2),c(1,1),col=8,lwd=7/12,lty=3)
points(res.brtr.r~pre.brtr,smg,col=col,pch=pch,bg=bg,cex=.5,lwd=.2)
points(res.brtr.r~pre.brtr,smg[smg$sp==2&smg$site=='dk',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)
points(res.brtr.r~pre.brtr,smg[smg$sp==2&smg$site=='fr',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)	
points(res.brtr.r~pre.brtr,smg[smg$sp==1&smg$site=='ok',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)	
mtext('b',line=-.8,adj=0.02,font=2,cex=8/12)
axis (1,seq(-2,2,by=.5),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-2,5,by=1),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Normalized residuals')),2,line=.7,font=1,cex=7/12)


# plot(res.brtr~pre.brtr,smg,ylim=c(-1,1),xlim=c(0,.8),xlab='Predicted BT using a 4-p logistic',ylab='Relative residuals',pch=pch,bg=bg,col=col)
# c("P. abies","E. grandis","P. taeda","P. sylvestris","P. menziesii")
## R-square
summary(lm(h.d.res.r~rsd.pre*as.factor(sp),smg))$r.s+(1-summary(lm(h.d.res.r~rsd.pre*as.factor(sp),smg))$r.s)*summary(lm(brtr~h.d.pre*as.factor(sp),smg))$r.s


################################################################################################
## Impact at stand scale
## models for NPPaw with/without branch turnover
################################################################################################
smg1<-smg[smg$sp==1,]
smg2<-smg[smg$sp==2,]

smg1$res.rc<-residuals(lm(rNPPbr~rc,smg1))
summary(lm(rNPPbr~rc*as.factor(sp),smg))
summary(lm(rNPPaw~rc*as.factor(sp),smg))

rNPPbr.rc<-lm(rNPPbr~rc+as.factor(sp),smg)
summary(rNPPbr.rc)
rNPPaw.rc<-lm(rNPPaw~rc+as.factor(sp),smg)
summary(rNPPaw.rc)

md(smg$brtr)
range(smg$brtr)
md(smg$rNPPbr)
range(smg$rNPPbr)

f.rNPPbr.rc<-function(rc,sp){
	x<-rc
	coe<-summary(rNPPbr.rc)$coe
	coe[1]+coe[3]*(sp-1)+coe[2]*(x)}
f.rNPPaw.rc<-function(rc,sp){
	x<-rc
	coe<-summary(rNPPaw.rc)$coe
	coe[1]+coe[3]*(sp-1)+(coe[2])*(x)}


Brtr.mod1<-data.frame(summary(nls(brtr~c+((d-c)/(1+exp(-a*(h.d-b)))),smg[smg$sp==1,],start=list(a=2.22,b=2.55,c=.09,d=.81)))$coe)
Brtr.mod1$mod<-'Brtr.mod1'
Brtr.res.mod11<-data.frame(summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))$coe)
Brtr.res.mod11$mod<-'Brtr.res.mod11'
Brtr.res.mod12<-data.frame(summary(nls(h.d.res.r~c+((d-c)/(1+exp(-a*(rsd-b)))),smg[smg$sp==1,],start=list(a=6.06,b=.7212,c=.5796,d=1.6384)))$coe)
Brtr.res.mod12$mod<-'Brtr.res.mod12'



Brtr.mod2<-data.frame(summary(lm(brtr~h.d,smg[smg$sp==2,]))$coe)
Brtr.mod2$mod<-'Brtr.mod2'
Brtr.res.mod2<-data.frame(summary(lm(h.d.res.r~std,smg[smg$sp==2,]))$coe)
Brtr.res.mod2$mod<-'Brtr.res.mod2'

rNPPbr.mod<-data.frame(summary(rNPPbr.rc)$coe)
rNPPbr.mod$mod<-'rNPPbr.mod'

rNPPaw.mod<-data.frame(summary(rNPPaw.rc)$coe)
rNPPaw.mod$mod<-'rNPPaw.mod'

parameters.scaled<-rbind(Brtr.mod1, Brtr.res.mod11, Brtr.res.mod12, Brtr.mod2, Brtr.res.mod2, rNPPbr.mod, rNPPaw.mod)

write.table(parameters.scaled,'/Users/hyli0001/Documents/wd/7_Branch_turnover/Figure/BT_Stand_scaled_parameters.txt',sep='\t',quote=FALSE,row.names=TRUE)
# # capture.output(parameters.scaled,file='/Users/hyli0001/Documents/wd/7_Branch_turnover/Figure/

### FIGURES
quartz(w=2,h=4.5)
par(lwd=.01)
par(mfrow=c(1,2))
par(mai=c(0,0,0,0))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,5),ylim=c(0,11))

lgr<-merge(sym,lnd[,c('site','st')],by='site')
lgdt<-lgr[order(lgr$ord),]
lgdt$strt<-paste(lgdt$st,lgdt$trt1,lgdt$trt2,sep='-')
lgd1<-lgdt[lgdt$sp==1,]
legend('topleft',lgd1$strt,pch=lgd1$pch,cex=0.4,col= lgd1$col,pt.bg=lgd1$bg,box.col=0,ncol=1,pt.lwd=.1,horiz=FALSE)

par(mai=c(0,0,0,0))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,5),ylim=c(0,11))
lgd2<-lgdt[lgdt$sp==2,]
legend('topleft',lgd2$strt,pch=lgd2$pch,cex=0.4,col= lgd2$col,pt.bg=lgd2$bg,box.col=0,ncol=1,pt.lwd=lgd2$lwd,horiz=FALSE)


# quartz(w=3.4,h=1.2)
# par(lwd=.01)
# par(mfrow=c(1,2))
# par(mai=c(0,.2,.5,0))
# plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,5.2),ylim=c(0,1.1))

# lnd1<-rbind(lnd[lnd$sp==1,][1:15,])
# mtext('Shade-intolerant',3,adj=.05,cex=0.5,line=-0.3)
# legend('topleft',x.intersp=.5,y.intersp=.85,lnd1$st,pch=21,cex=6/12,col=lnd1$col,pt.bg=lnd1$bg,box.col=0,ncol=5,pt.lwd=.1)

# par(mai=c(0,0,.5,.2))
# plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,5.2),ylim=c(0,1.1))
# lnd2<-lnd[lnd$sp==2,]
# mtext('Shade-tolerant',3,adj=.05,cex=0.5,line=-0.3)
# legend('topleft',x.intersp=.5,y.intersp=.85,lnd2$st,pch=21,cex=6/12,col=lnd2$col,pt.bg=lnd2$bg,box.col=0,ncol=3,pt.lwd=.75)


quartz(w=3.42,h=3.6)
par(mfrow=c(2,2))
par(lwd=.3)
par(mai=c(.2,.4,.4,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,5),ylim=c(0,1))
points(brtr~h.d,smg[smg$sp==1,],col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
mtext('a',line=-.8,adj=0.02,font=2,cex=8/12)
mtext('Shade-intolerant species',line=-.1,adj=1,font=1,cex=7/12)
summary(lm(h.d.pre~brtr,smg[smg$sp==1,]))
mtext(expression(paste('r'^2,'=0.91; p<.001',sep='')),1,line=-1,adj=1,font=1,cex=7/12)

axis (1,seq(-2,9,by=2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-1,2,by=.5),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Annual turnover rate of branch biomass (y'^-1,')')),2,line=.6,font=1,cex=7/12,at=-.2)
# mtext(expression(paste('(kg T'[BR],' kg'^-1,' W'[BR],' yr'^-1,')')),2,line=.5,font=1,cex=6/12)
mtext(expression(paste('Height increment (m y'^-1,')')),1,line=.5,font=1,cex=7/12)
unique(df$site[df$site!='hb'])
for(i in 1){
	df<-smg[smg$sp==i,]
	site<-unique(df$site[df$site!='hb'])
	site<-site[13:1]
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(brtr+h.d~bg+pch+col,df2,FUN=me)
			lines(rep(df3$h.d.m,2),c(-1,1)*df3$brtr.se+df3$brtr.m)
			lines(c(-1,1)*df3$h.d.se+df3$h.d.m,rep(df3$brtr.m,2))
			points(brtr.m~h.d.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=0.2)
			}}}
		points(brtr~h.d,smg[smg$site=='hb',],bg=bg,pch=pch,col=col,cex=7/12,lwd=0.2)
	curve(h.d(x,1),add=TRUE,xlim=range(smg$h.d[smg$sp==1],na.rm=TRUE),lwd=0.75)
lnd1<-rbind(lnd[lnd$sp==1,][1:9,])
legend(-.3,.95,x.intersp=.5,y.intersp=.85,lnd1$st,pch=21,cex=7/12,col=lnd1$col,pt.bg=lnd1$bg,box.col=0,ncol=3,pt.lwd=.1)
lnd2<-rbind(lnd[lnd$sp==1,][10:15,])
legend(-.3,.73,x.intersp=.5,y.intersp=.85,lnd2$st,pch=21,cex=7/12,col=lnd2$col,pt.bg=lnd2$bg,box.col=0,ncol=2,pt.lwd=.1)

par(mai=c(.2,.35,.4,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(-0.1,2.5),xlim=c(0,1.45))
points(h.d.res.r~rsd,smg[smg$sp==1,],col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
mtext('b',line=-.8,adj=0.02,font=2,cex=8/12)
summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))
mtext(expression(paste('r'^2,'=0.51; p<.001')),1,line=-1,adj=1,font=1,cex=7/12)
axis (1,seq(-2,8,by=.5),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-2,8,by=1),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Relative residuals (observed / predicted)')),2,line=.6,font=1,cex=7/12,at=-.5)
mtext(expression(paste('Relative stand density')),1,line=.5,font=1,cex=7/12)
for(i in 1){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df2$ln.std<-log(df2$std)
			df3<-summaryBy(h.d.res.r+ln.std+std+rsd~bg+pch+col,df2,FUN=me)
			lines(rep(df3$rsd.m,2),c(-1,1)*df3$h.d.res.r.se+df3$h.d.res.r.m)
			lines(c(-1,1)*df3$rsd.se+df3$rsd.m,rep(df3$h.d.res.r.m,2))
			points(h.d.res.r.m~rsd.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=0.2)
			}}}
curve(h.d.std.r(NA,x,1,2),xlim=range(smg$rsd[smg$sp==1],na.rm=TRUE),add=TRUE,lwd=0.75)

par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,0.9),ylim=c(-.01,.09))
mtext('c',line=-.8,adj=0.02,font=2,cex=8/12)
mtext('Shade-tolerant species',line=-.1,adj=1,font=1,cex=7/12)
summary(lm(brtr~h.d,smg[smg$sp==2,]))
mtext(expression(paste('r'^2,'=0.43; p<.001',sep='')),1,line=-1,adj=1,font=1,cex=7/12)
points(brtr~h.d,smg[smg$sp==2,],col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
axis (1,seq(-2,9,by=0.2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-0.4,1,by=0.04),,tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Height increment (m y'^-1,')')),1,line=.5,font=1,cex=7/12)

for(i in 2){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(brtr+h.d~bg+pch+col,df2,FUN=me)
			lines(rep(df3$h.d.m,2),c(-1,1)*df3$brtr.se+df3$brtr.m)
			lines(c(-1,1)*df3$h.d.se+df3$h.d.m,rep(df3$brtr.m,2))
			points(brtr.m~h.d.m,df3,bg=bg,pch=pch,col=col,cex=5/12,lwd=0.3)
			points(brtr.m~h.d.m,df3,bg=0,pch=pch,col=1,cex=7/12,lwd=0.3)
			}}}
curve(h.d(x,2),add=TRUE,xlim=range(smg$h.d[smg$sp==2],na.rm=TRUE),lwd=0.65)
lnd3<-lnd[lnd$sp==2,]
legend(-.05,0.085,x.intersp=.5,y.intersp=.85,lnd3$st,pch=21,cex=7/12,col=lnd3$col,pt.bg=lnd3$bg,box.col=0,ncol=3,pt.lwd=.75)


par(mai=c(.4,.35,.2,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,2.5),xlim=c(250,3400))
mtext('d',line=-.8,adj=0.02,font=2,cex=8/12)
summary(lm(h.d.res.r~std,smg[smg$sp==2,]))
mtext(expression(paste('r'^2,'=0.41; p<.001',sep='')),1,line=-1,adj=1,font=1,cex=7/12)
points(h.d.res.r~std,smg[smg$sp==2,],col=col,pch=pch,bg='white',cex=0.4,lwd=.05)
axis (1,seq(0,6000,by=1000),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-1,5,by=1),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Stand density (stems ha'^-1,')')),1,line=.5,font=1,cex=7/12)

for(i in 2){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(h.d.res.r+sdi+std~bg+pch+col,df2,FUN=me)
			lines(rep(df3$std.m,2),c(-1,1)*df3$h.d.res.r.se+df3$h.d.res.r.m)
			lines(c(-1,1)*df3$std.se+df3$std.m,rep(df3$h.d.res.r.m,2))
			points(h.d.res.r.m~std.m,df3,bg=bg,pch=pch,col=col,cex=5/12,lwd=0.3)
			points(h.d.res.r.m~std.m,df3,bg=0,pch=pch,col=1,cex=7/12,lwd=0.3)
			}}}
curve(h.d.std.r(x,NA,2,2),xlim=range(smg$std[smg$sp==2],na.rm=TRUE),add=TRUE,lwd=0.65)


quartz(w=2,h=4.5)
par(lwd=.01)
par(mfrow=c(1,2))
par(mai=c(0,0,0,0))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,5),ylim=c(0,11))

lgr<-merge(sym,lnd[,c('site','st')],by='site')
lgdt<-lgr[order(lgr$ord),]
lgdt$strt<-paste(lgdt$st,lgdt$trt1,lgdt$trt2,sep='-')
lgd1<-lgdt[lgdt$sp==1,]
legend('topleft',lgd1$strt,pch=lgd1$pch,cex=0.4,col= lgd1$col,pt.bg=lgd1$bg,box.col=0,ncol=1,pt.lwd=.1,horiz=FALSE)

par(mai=c(0,0,0,0))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,5),ylim=c(0,11))
lgd2<-lgdt[lgdt$sp==2,]
legend('topleft',lgd2$strt,pch=lgd2$pch,cex=0.4,col= lgd2$col,pt.bg=lgd2$bg,box.col=0,ncol=1,pt.lwd=lgd2$lwd,horiz=FALSE)


quartz(w=3.42,h=1.8)
par(mfrow=c(1,2))
par(lwd=.3)
par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(.15,1),ylim=c(0,1.25))
points(rNPPbr~rc,smg,col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
mtext('a',line=-.8,adj=0.02,font=2,cex=8/12)
summary(rNPPbr.rc)
mtext(expression(paste('r'^2,'=0.56; p<.001')),1,line=-.9,adj=.05,font=1,cex=7/12)
mtext(expression(paste('')),line=-1.3,adj=.95,font=1,cex=5/12)
axis (1,seq(-2,2,by=.2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-2,2,by=.5),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' : NPP'[B])),2,line=.6,font=1,cex=7/12)
mtext(expression(paste('Live crown ratio')),1,line=.5,font=1,at=1.1,cex=7/12)

for(i in 1:2){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(rNPPbr+rc~bg+pch+col,df2,FUN=me)
			lines(rep(df3$rc.m,2),c(-1,1)*df3$rNPPbr.se+df3$rNPPbr.m)
			lines(c(-1,1)*df3$rc.se+df3$rc.m,rep(df3$rNPPbr.m,2))
			points(rNPPbr.m~rc.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=i*0.1+0.1)
			if(i==2){
			points(rNPPbr.m~rc.m,df3,bg=0,pch=pch,col=1,cex=7.5/12,lwd=0.1)}
			}}
curve(f.rNPPbr.rc(x,i),add=TRUE,xlim=range(smg$rc[smg$sp==i],na.rm=TRUE),lwd=0.75)}

par(mai=c(.4,.35,.2,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(.1,1),ylim=c(0,.3))
points(rNPPaw~rc,smg,col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
mtext('b',line=-.8,adj=0.02,font=2,cex=8/12)
summary(rNPPaw.rc)
mtext(expression(paste('r'^2,'=0.40; p<.001')),1,line=-.9,adj=.05,font=1,cex=7/12)
mtext(expression(paste('')),line=-1.3,adj=.95,font=1,cex=5/12)
axis (1,seq(-2,2,by=.2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-2,2,by=.1),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' : NPP'[Wa])),2,line=.6,font=1,cex=7/12)

for(i in 1:2){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(rNPPaw+rc~bg+pch+col,df2,FUN=me)
			lines(rep(df3$rc.m,2),c(-1,1)*df3$rNPPaw.se+df3$rNPPaw.m)
			lines(c(-1,1)*df3$rc.se+df3$rc.m,rep(df3$rNPPaw.m,2))
			points(rNPPaw.m~rc.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=i*0.1+0.1)
			if(i==2){
			points(rNPPaw.m~rc.m,df3,bg=0,pch=pch,col=1,cex=7.5/12,lwd=0.1)}
			}}
curve(f.rNPPaw.rc(x,i),add=TRUE,xlim=range(smg$rc[smg$sp==i],na.rm=TRUE),lty=1,lwd=0.75)}


## Gamma distribution : shape=mean^2/var; rate=mean/var; scale=var/mean; Dispersal=var/mean^2

quartz(w=4.5,h=2.55)
par(mfrow=c(1,4))
par(lwd=.1,col=0)
par(mai=c(.35,0,.2,.0))
plot(NA,xlab="",ylab="",yaxt="n",xaxt="n",xlim=c(0,10),ylim=c(-.2,9.15))
# ,ylim=c(.5,8.35)
x.label<-unique(smg[,c('spe','rev.ord')])
x.label<-x.label[order(-x.label$rev.ord),]
x.label$ful.nm<-c('Pinus sylvestris','Pinus taeda','Eucalyptus grandis','Pseudotsuga menziesii','Populus tremula × tremuloides','Betula pendula','Bruguiera gymnorrhiza','Picea abies','Broad-leaved deciduous')

par(lwd=.1,col=1)
mtext(x.label$ful.nm[1:8],las=2,2,at=(x.label$rev.ord[1:8])-0.1,line=-8.5,font=3,cex=6/12)
mtext(x.label$ful.nm[9],las=2,2,at=(x.label$rev.ord[9])-0.1,line=-8.5,font=c(1),cex=6/12)
mtext('Average',las=2,2,at=-0.25,line=-8,font=c(1),cex=6/12)

for (i in c(6:9)){
	lnb<-rbind(lnd[1,],lnd[lnd$rev.ord%in%i,])
	lnb[1,]<-NA
legend(yjust=1.25,3,i,x.intersp=.4,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}
for (i in c(3:5)){
	lnb<-rbind(lnd[1:3,],lnd[lnd$rev.ord%in%i,])
lnb[1:3,]<-NA
legend(yjust=1.25,3,i,x.intersp=.4,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}

for (i in c(1)){
	lnb<-rbind(lnd[1:3,],lnd[lnd$rev.ord%in%i,])
lnb[1:3,]<-NA
legend(yjust=1.25,3,i,x.intersp=.4,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=1)
}

for (i in c(2)){
lnb<-lnd[lnd$rev.ord%in%i,][1:7,]
legend(yjust=1.2,-.3,i,x.intersp=.4,y.intersp=.85,lnb$st,cex=7/12,pch=21,col=lnb$col,pt.bg=lnb$bg,box.col=0,pt.lwd=ifelse(lnb$sp==1,.2,.7),horiz=TRUE,text.width=.8)
}


par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,0.9),ylim=c(-.2,9.25))
points((jitter(rev.ord,.5)-0.5)~brtr,smg,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(smg$sp==1,bg,col))
mtext('a',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-4,4,by=0.5),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('∆'[B],' (y'^-1,')')),1,line=.9,font=1,cex=7/12)

for (i in 1:2){
	df<-smg[smg$sp==i,]
	for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$brtr[dff$brtr>=0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,2,0.001),shape,rate))*1.85
	y.mn<-unique(dff$rev.ord)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
	# print(unique(dff[,c('site','sp')]))
	# print(gamma_test(x))
}}

	for (i in 1:nrow(smm)){
		lines(c(-1,1)*smm$brtr.sd[i]+smm$brtr.m[i],rep(smm$rev.ord[i],2)-.5)}
	points((rev.ord-.5)~brtr.m,smm,bg=bg,pch=21,col=col,cex=7/12,lwd=ifelse(smm$sp==1,.2,.7))


xv<-smg$brtr[smg$brtr>=0]
shape<-mean(xv)^2/var(xv)
rate<-mean(xv)/var(xv)
t.max<-max(dnorm(seq(0.01,1,0.001),mean(xv),var(xv)))*0.85
points(mean(xv),-.3,pch=21,bg=1,cex=8/12)
lines(c(-1,1)*sqrt(var(xv))+mean(xv),rep(-.3,2))
# curve(dgamma(x,shape,rate)/t.max-.3,col=1,lwd=.7,add=TRUE,xlim=c(range(xv)[1]*3,range(xv)[2]))


par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,1.1),ylim=c(-.2,9.15))
points((jitter(rev.ord,.5)-0.5)~rNPPbr,smg,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(smg$sp==1,bg,col))
mtext('b',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-4,4,by=0.5),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' : NPP'[B])),1,line=.9,font=1,cex=7/12)

for (i in 1:2){
df<-smg[smg$sp==i,]
for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$rNPPbr[dff$rNPPbr>=0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,2,0.001),shape,rate))*1.85
	y.mn<-unique(dff$rev.ord)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
	# print(unique(dff[,c('site','sp')]))
	# print(gamma_test(x))
	}}

	for (i in 1:nrow(smm)){
		lines(c(-1,1)*smm$rNPPbr.sd[i]+smm$rNPPbr.m[i],rep(smm$rev.ord[i],2)-.5)}
	points((rev.ord-.5)~rNPPbr.m,smm,bg=bg,pch=21,col=col,cex=7/12,lwd=ifelse(smm$sp==1,.2,.7))

xv<-smg$rNPPbr[dff$rNPPbr>=0]
shape<-mean(xv)^2/var(xv)
rate<-mean(xv)/var(xv)
t.max<-max(dnorm(seq(0.01,1,0.001),mean(xv),var(xv)))*0.85
points(mean(xv),-.3,pch=21,bg=1,cex=8/12)
lines(c(-1,1)*sqrt(var(xv))+mean(xv),rep(-.3,2))
# curve(dgamma(x,shape,rate)/t.max-.3,col=1,lwd=.7,add=TRUE,xlim=range(xv))

par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,.31),ylim=c(-.2,9.15))
points(jitter(rev.ord,.5)-.5~rNPPaw,smg,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(smg$sp==1,bg,col))
mtext('c',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-1,1,by=.1),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' : NPP'[Wa])),1,line=.9,font=1,cex=7/12)
mtext('Scaled plots',line=-.1,adj=1,font=1,cex=7/12)


for (i in 1:2){
df<-smg[smg$sp==i,]
for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$rNPPaw[dff$brtr>0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,2,0.001),shape,rate))*1.85
	y.mn<-unique(dff$rev.ord)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
	# print(unique(dff[,c('site','sp')]))
	# print(gamma_test(x))
	}}

	for (i in 1:nrow(smm)){
		lines(c(-1,1)*smm$rNPPaw.sd[i]+smm$rNPPaw.m[i],rep(smm$rev.ord[i],2)-.5)}
	points((rev.ord-.5)~rNPPaw.m,smm,bg=bg,pch=21,col=col,cex=7/12,lwd=ifelse(smm$sp==1,.2,.7))

xv<-smg$rNPPaw[dff$rNPPaw>=0]
shape<-mean(xv)^2/var(xv)
rate<-mean(xv)/var(xv)
t.max<-max(dnorm(seq(0.01,1,0.001),mean(xv),var(xv)))*0.85
points(mean(xv),-.3,pch=21,bg=1,cex=8/12)
lines(c(-1,1)*sqrt(var(xv))+mean(xv),rep(-.3,2))
# curve(dgamma(x,shape,rate)/t.max-.3,col=1,lwd=.7,add=TRUE,xlim=range(xv))


####
########### Litter fall #######

dklrr<-read.table('/Users/hyli0001/Documents/wd/7_Branch_turnover/Raw_data/DukeFACE/Rinput_Leaf litter data.csv',head=TRUE,sep=',')
dklrr$date<-as.Date(dklrr$Collection_date,'%m/%d/%y')
dkltr<-dklrr[dklrr$Year%in%2005:2013,]
dkltr$Pine_branch_woodonly<-NULL
dkltr$Hardwood_branch_woodonly<-NULL
dkltr$All_barkonly<-NULL


dkltr$mon<-ifelse(dkltr$Month=='JAN',1,ifelse(dkltr$Month=='FEB',2,ifelse(dkltr$Month=='MAR',3,ifelse(dkltr$Month=='APR',4,ifelse(dkltr$Month=='MAY',5,ifelse(dkltr$Month=='JUN',6,ifelse(dkltr$Month=='JUL',7,ifelse(dkltr$Month=='AUG',8,ifelse(dkltr$Month=='SEP',9,ifelse(dkltr$Month=='OCT',10,ifelse(dkltr$Month=='NOV',11,ifelse(dkltr$Month=='DEC',12,NA))))))))))))

dkltr$branch<-dkltr$Pine_branch+dkltr$Hardwood_branch
dlym<-summaryBy(branch+Pine_branch+Hardwood_branch+Pine_leaf_dead+Hardwood_leaf+Pine_cone+Pine_seed+Other_seed+Fine_debris~Plot+Quadrant+Year+mon,dkltr,FUN=sum,keep.names=TRUE)

# dlym[dlym$Pine_branch>500,]
### Giving a threshold of ~340 g m2 ±3SD
otl<-mean(dlym$branch)+sd(dlym$branch)*3

dlym$lf_br<-ifelse(dlym$branch<otl,dlym$branch,NA)
dlym$branch.avg<-ave(dlym$lf_br,dlym$Plot,dlym$Quadrant,dlym$mon,FUN=function(x){mean(x,na.rm=TRUE)})
dlym$lf_br<-ifelse(is.na(dlym$lf_br),dlym$branch.avg,dlym$lf_br)
dlymm<-summaryBy(lf_br~Year+mon,dlym,FUN=mean,keep.names=TRUE)

# plot(Pine_branch~mon,dlym[dlym$Year==2006,],ylim=c(0,1500))
# plot(Pine_branch~mon,dlym[dlym$Year==2007,],ylim=c(0,1500))
# plot(Pine_branch~mon,dlym[dlym$Year==2008,],ylim=c(0,1500))
# plot(Pine_branch~mon,dlym[dlym$Year==2009,],ylim=c(0,1500))
# plot(Pine_branch~mon,dlym[dlym$Year==2010,],ylim=c(0,1500))

dklyr<-summaryBy(lf_br~Plot+Quadrant+Year,dlym,FUN=sum,keep.names=TRUE)
duke.plot.info<-unique(read.table('/Users/hyli0001/Documents/wd/7_Branch_turnover/Raw_data/DukeFACE/TreeInfo_WithoutGapfilledHeights.csv',sep=',',head=TRUE)[,3:6])
colnames(duke.plot.info)[2]<-c('Quadrant')

dkl<-merge(dklyr,duke.plot.info,by=c('Plot','Quadrant'))
dkl$Quadrant<-NULL
dkl$trt1<-ifelse(dkl$CO2=='amb','A','E')
dkl$CO2<-NULL
dkl$trt2<-ifelse(dkl$N=='cont','C','F')
dkl$N<-NULL
colnames(dkl)<-c('plot','year','lf_br','trt1','trt2')
dklb<-summaryBy(lf_br~plot+year+trt1+trt2, dkl ,FUN=mean,keep.names=TRUE)
dkp<-summaryBy(NPPbr0+brt+dBdb~year+plot+trt1+trt2,btm[btm$site=='dk',],FUN=sum,keep.names=TRUE)
dkp$bg<-ifelse(dkp$trt1=='E',2,4)
dkp$pch<-ifelse(dkp$trt2=='F',24,21)
dkp$tt<-paste(dkp$trt1,dkp$trt2,sep='')
dkp$NPPbr<-dkp$NPPbr0+dkp$brt
dkp$brt.adj<-dkp$brt/0.840
dkp$NPPbr.adj<-dkp$NPPbr0+dkp$brt.adj

dkpl<-merge(dkp,dklb[dklb$year%in%2006:2011&dklb$plot%in%1:8,],by=c('plot','year','trt1','trt2'))
dkpl$tt<-paste(dkpl$trt1,dkpl$trt2,sep='')
dkpl$x.lab<-ifelse(dkpl$tt=='AC',1,ifelse(dkpl$tt=='AF',2,ifelse(dkpl$tt=='EC',3,ifelse(dkpl$tt=='EF',4,NA))))
dkpl$mcfr<-1-((dkpl$lf_br+dkpl$dBdb)/dkpl$brt)
dkpl$mcfr.adj<-1-((dkpl$lf_br+dkpl$dBdb)/dkpl$brt.adj)
dkpl$mcfa<-(dkpl$brt-dkpl$lf_br-dkpl$dBdb)
dkpl$mcfa.adj<-(dkpl$brt.adj-dkpl$lf_br-dkpl$dBdb)
dkpt<-summaryBy(mcfr+mcfa+mcfr.adj+mcfa.adj+Bbr+brt+lf_br+dBdb~trt1+trt2+tt+year+bg+x.lab+pch,dkpl,FUN=me)
dkpm<-summaryBy(NPPbr+brt+lf_br+dBdb~year,dkpl[dkpl$tt=='AC',],FUN=me)

quartz(w=3.42,h=3.6)
par(mfrow=c(2,2))
par(lwd=.3)
par(mai=c(.1,.4,.5,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,950),xlim=c(2005.5,2010.5))
points(NPPbr~year,dkpl[dkpl$tt=='AC',],bg=8,pch=21,cex=0.7,lwd=0.05)
mtext('a',line=-.8,adj=0.05,font=2,cex=8/12)
mtext('Pinus taeda',line=-.1,adj=0,font=3,cex=7/12)
mtext('& broad-leaved deciduous',line=-.1,adj=3.4,font=1,cex=7/12)

for (i in 1:5){
	lines(rep(dkpm$year[i],2),dkpm$NPPbr.m[i]+c(-1,1)*dkpm$NPPbr.se[i])
	}
	points(NPPbr.m~year,dkpm,bg='red4',pch=21,cex=0.75)
coe1<-summary(lm(NPPbr~year,dkpl[dkpl$tt=='AC',]))$coe
curve(coe1[1]+coe1[2]*x,add=TRUE,xlim=c(2006,2010))
mtext('p = 0.006',line=-.8,adj=0.9,font=1,cex=7/12)
axis (1,seq(2004,2020,by=2),tck=.02,label=FALSE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-600,1200,by=300),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('NPP'[B],' (g C m'^-2,' y'^-1,')')),2,line=.6,font=1,cex=7/12)

par(mai=c(.1,.35,.5,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,750),xlim=c(2005.5,2010.5))
points(brt~year,dkpl[dkpl$tt=='AC',],bg=8,pch=21,cex=0.7,lwd=0.05)
mtext('b',line=-.8,adj=0.05,font=2,cex=8/12)
mtext('Duke FACE',line=-.1,adj=1,font=1,cex=7/12)
for (i in 1:5){
	lines(rep(dkpm$year[i],2),dkpm$brt.m[i]+c(-1,1)*dkpm$brt.se[i])
	}
	points(brt.m~year,dkpm,bg='red4',pch=21,cex=0.75)
coe2<-summary(lm(brt~(year),dkpl[dkpl$tt=='AC',]))$coe
curve(coe2[1]+coe2[2]*x,add=TRUE,xlim=c(2006,2010))
mtext('p = 0.008',line=-.8,adj=0.9,font=1,cex=7/12)

axis (1,seq(2004,2020,by=2),tck=.02,label=FALSE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-600,1200,by=200),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' (g C m'^-2,' y'^-1,')')),2,line=.6,font=1,cex=7/12)

par(mai=c(.5,.4,.1,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,450),xlim=c(2005.5,2010.5))
points(lf_br~year,dkpl[dkpl$tt=='AC',],bg=8,pch=21,cex=0.7,lwd=0.05)
mtext('c',line=-.8,adj=0.05,font=2,cex=8/12)
for (i in 1:5){
	lines(rep(dkpm$year[i],2),dkpm$lf_br.m[i]+c(-1,1)*dkpm$lf_br.se[i])
	}
	points(lf_br.m~year,dkpm,bg='red4',pch=21,cex=0.75)
summary(lm(lf_br~year,dkpl[dkpl$tt=='AC',]))
curve(mean(dkpl[dkpl$tt=='AC','lf_br'])+0*x,add=TRUE,xlim=c(2006,2010),lty=3)
mtext('p = 0.242',line=-.8,adj=0.9,font=1,cex=7/12)

axis (1,seq(2004,2020,by=2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-600,1200,by=200),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('F'[B],' (g C m'^-2,' y'^-1,')')),2,line=.5,font=1,cex=7/12)
mtext(expression(paste('Year')),1,line=.5,font=1,cex=7/12)

par(mai=c(.5,.35,.1,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,55),xlim=c(2005.5,2010.5))
points(dBdb~year,dkpl[dkpl$tt=='AC',],bg=8,pch=21,cex=0.7,lwd=0.05)
mtext('d',line=-.8,adj=0.05,font=2,cex=8/12)
for (i in 1:5){
	lines(rep(dkpm$year[i],2),dkpm$dBdb.m[i]+c(-1,1)*dkpm$dBdb.se[i])
	}
	points(dBdb.m~year,dkpm,bg='red4',pch=21,cex=0.75)
summary(lm(dBdb~year,dkpl[dkpl$tt=='AC',]))
curve(mean(dkpl[dkpl$tt=='AC','dBdb'])+0*x,add=TRUE,xlim=c(2006,2010),lty=3)
mtext('p = 0.923',line=-.8,adj=0.9,font=1,cex=7/12)

axis (1,seq(2004,2020,by=2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-600,1200,by=20),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('∆W'[DB],' (g C m'^-2,' y'^-1,')')),2,line=.5,font=1,cex=7/12)
mtext(expression(paste('Year')),1,line=.5,font=1,cex=7/12)

dkp$bg<-ifelse(dkp$trt1=='E',2,4)
dkp$pch<-ifelse(dkp$trt2=='F',24,21)

dbd.p<-read.table('/Users/hyli0001/Documents/wd/7_Branch_turnover/Raw_data/DukeFACE/Fig2S_pine.csv',head=TRUE,sep=',')
dbd.p$bg<-ifelse(dbd.p$Treatment%in%c('EC','EF'),2,4)
dbd.p$pch<-ifelse(dbd.p$Treatment%in%c('AF','EF'),24,21)

dbd.h<-read.table('/Users/hyli0001/Documents/wd/7_Branch_turnover/Raw_data/DukeFACE/Fig2S_hw.csv',head=TRUE,sep=',')
dbd.h$bg<-ifelse(dbd.h$Treatment%in%c('EC','EF'),2,4)
dbd.h$pch<-ifelse(dbd.h$Treatment%in%c('AF','EF'),24,21)

quartz(w=3.42,h=3.6)
par(mfrow=c(2,2))
par(lwd=.3)
par(mai=c(.1,.4,.5,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(1.9,3.8),ylim=c(5,11.5))
mtext('a',line=-.8,adj=0.05,font=2,cex=8/12)
mtext('Pinus taeda',line=-.1,adj=0,font=3,cex=7/12)
mtext(expression(paste('r'^2,'=0.76; p<.001',sep='')),line=-.95,adj=.95,font=1,cex=7/12)
mtext(expression(paste('y = -1.97 + 3.29x',sep='')),line=-1.5,adj=.95,font=1,cex=7/12)
points(pine_Biomass_log~pine_DBH_log,dbd.p,col=1,pch=pch,bg=bg,cex=0.4,lwd=0.3)
curve(-1.97+3.29*x,add=TRUE,xlim=range(dbd.p$pine_DBH_log))
axis (1,seq(-2,9,by=.5),tck=.02,label=FALSE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-2,12,by=2),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('log (W'[DB],')')),2,line=.6,font=1,cex=7/12)
legend('bottomright',xjust=0.1,c('AC','AF','EC','EF'),col=1,pt.bg=c(4,4,2,2),pch=c(21,24,21,24),cex=7/12,box.lwd=0,text.font=1)

par(mai=c(.1,.35,.5,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,3.1),ylim=c(-2,11))
mtext('b',line=-.8,adj=0.05,font=2,cex=8/12)
mtext('Broad-leaved deciduous',line=-.1,adj=1,font=1,cex=7/12)
mtext(expression(paste('r'^2,'=0.34; p<.001',sep='')),line=-.95,adj=.95,font=1,cex=7/12)
mtext(expression(paste('y = 1.33 + 1.92x',sep='')),line=-1.5,adj=.95,font=1,cex=7/12)
points(hw_Biomass_log~hw_DBH_log,dbd.h,col=1,pch=pch,bg=bg,cex=0.4,lwd=.3)
curve(1.33+1.92*x,add=TRUE,xlim=range(dbd.h$hw_DBH_log))
axis (1,seq(-2,9,by=1),tck=.02,label=FALSE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-4,12,by=4),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('log (W'[DB],')')),2,line=.6,font=1,cex=7/12)

par(mai=c(.5,.4,.1,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(1.9,3.8),ylim=c(-4,4))
mtext('c',line=-.8,adj=0.05,font=2,cex=8/12)
points(pine_stdresid~pine_DBH_log,dbd.p,col=1,pch=pch,bg=bg,cex=0.4,lwd=.3)
curve(0*x,add=TRUE,lty=3)
axis (1,seq(-2,9,by=.5),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-6,6,by=2),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Standardized residuals')),2,line=.6,font=1,cex=7/12)
mtext(expression(paste('log (diamater at 1.3 m)')),1,line=.6,font=1,cex=7/12)

par(mai=c(.5,.35,.1,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,3.1),ylim=c(-4,4))
mtext('d',line=-.8,adj=0.05,font=2,cex=8/12)
points(hw_stdresid~hw_DBH_log,dbd.h,col=1,pch=pch,bg=bg,cex=0.4,lwd=.3)
curve(0*x,add=TRUE,lty=3)
axis (1,seq(-2,9,by=1),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-6,6,by=2),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Standardized residuals')),2,line=.6,font=1,cex=7/12)
mtext(expression(paste('log (diamater at 1.3 m)')),1,line=.6,font=1,cex=7/12)

mean(dkpl[dkpl$tt=='AC',]$dBdb/dkpl[dkpl$tt=='AC',]$lf_br)
mean(dkpl[dkpl$tt=='AC',]$dBdb/dkpl[dkpl$tt=='AC',]$brt)
mean(dkpl[dkpl$tt=='AC',]$lf_br/dkpl[dkpl$tt=='AC',]$brt)
	

dkps<-summaryBy(NPPbr0+NPPbr+NPPbr.adj+dBdb+brt+brt.adj+lf_br~trt1+trt2+plot+tt+bg+x.lab+pch,dkpl[dkpl$tt=='AC',],FUN=function(x){mean(x,na.rm=T)/2},keep.names=TRUE)
dkps$mcfa<-(dkps$brt-(dkps$lf_br+dkps$dBdb))
dkps$mcfa.adj<-(dkps$brt.adj-(dkps$lf_br+dkps$dBdb))
dkps$mcfr<-1-((dkps$lf_br+dkps$dBdb)/dkps$brt)
dkps$mcfr.adj<-1-((dkps$lf_br+dkps$dBdb)/dkps$brt.adj)


dkpsm<-summaryBy(NPPbr+brt+lf_br+dBdb+mcfa+mcfr+NPPbr0~site,dkps,FUN=me)
dkpsm.adj<-summaryBy(NPPbr.adj+brt.adj+lf_br+dBdb+mcfa.adj+mcfr.adj~site,dkps,FUN=me)
me(dkps$brt.adj/dkps$NPPbr.adj)
quartz(w=3.42,h=1.8)
par(mfrow=c(1,1))
par(lwd=.3)
par(mai=c(.4,.4,.2,.4))
plot(NA,xlim=c(0,11.5),ylim=c(9,300),xaxt='n',xlab='',ylab='', yaxt='n')
mtext('Branch carbon dynamics',line=-.1,adj=0.005,font=1,cex=7/12)
mtext('Duke FACE, 2006-2010',line=-.1,adj=1,font=1,cex=7/12)
df.adj<-dkpsm.adj
lines(c(1,4.25),rep(dkpsm[1],2),lty=3,lwd=1)
lines(c(2.5,4.25),rep(dkpsm[3],2),lty=3,lwd=1)
text(5.25,(dkpsm[1]+dkpsm[3])/2,'Increment',adj=.5,cex=7/12)

for(i in seq(1,9,2)){
	polygon(c(-.5,-.5,.5,.5)+i,c(-50,df.adj[1,i],df.adj[1,i],-50),col= 'darkorange1',border=1,lwd=.75)
	lines(rep(i,2), df.adj[1,i]+c(-1,1)* df.adj[1,i+1])}

for(i in seq(1,9,2)){
	polygon(c(-.5,-.5,.5,.5)+i,c(-50,dkpsm[1,i],dkpsm[1,i],-50),col='forestgreen',border=1,lwd=.75)
	lines(rep(i,2),dkpsm[1,i]+c(-1,1)*dkpsm[1,i+1])}

axis (1,seq(0,000,by=2000),tck=.02,label=c(expression(paste('NPP')),expression(paste('Turnover')),expression(paste('F')),expression(paste('∆W'[dead])),'Mass loss'),at=1:5*2-1,mgp=c(0,-.5,0),cex.axis=7/12,lwd=0.3,line=.35)

axis (2,seq(-200,800,by=100),tck=.02,label=TRUE,mgp=c(0,-.1,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Branch compartment')),1,line=.6,font=1,cex=7/12)
mtext(expression(paste('Carbon flux (g C m'^-2,' y'^-1,')')),2,line=.6,font=1,cex=7/12)

adj.tb<-dkpsm.adj$mcfa.adj.m-dkpsm$mcfa.m
ht<-275
polygon(c(-.3,-.3,.3,.3)+5.5,c(ht,ht+adj.tb,ht+adj.tb,ht),col= 'darkorange1',border=1,lwd=.75)
text(6,275,expression(paste('Within-canopy turnover')),adj=c(0,0),cex=7/12)

text(7.6,dkpsm$mcfa.m+60,expression(paste('Mass loss / Turnover')),adj=-.1,cex=6/12,font=2)
text(7.9,dkpsm$mcfa.m,'27 ± 16%',adj=-1,cex=6/12)
text(7.9,dkpsm.adj$mcfa.adj.m,'38 ± 13%',adj=-1,cex=6/12)


###### High & low ####
dkps.h<-summaryBy(NPPbr0+NPPbr+NPPbr.adj+dBdb+brt+brt.adj+lf_br~trt1+trt2+plot+tt+bg+x.lab+pch,dkpl[dkpl$tt=='AC'&dkpl$year%in%2006:2007,],FUN=function(x){mean(x,na.rm=T)/2},keep.names=TRUE)
dkps.h$mcfa<-(dkps.h$brt-(dkps.h$lf_br+dkps.h$dBdb))
dkps.h$mcfa.adj<-(dkps.h$brt.adj-(dkps.h$lf_br+dkps.h$dBdb))
dkps.h$mcfr<-1-((dkps.h$lf_br+dkps.h$dBdb)/dkps.h$brt)
dkps.h$mcfr.adj<-1-((dkps.h$lf_br+dkps.h$dBdb)/dkps.h$brt.adj)

dkpsm.h<-summaryBy(NPPbr+brt+lf_br+dBdb+mcfa+mcfr+NPPbr0~site,dkps.h,FUN=me)
dkpsm.adj.h<-summaryBy(NPPbr.adj+brt.adj+lf_br+dBdb+mcfa.adj+mcfr.adj~site,dkps.h,FUN=me)

dkps.l<-summaryBy(NPPbr0+NPPbr+NPPbr.adj+dBdb+brt+brt.adj+lf_br~trt1+trt2+plot+tt+bg+x.lab+pch,dkpl[dkpl$tt=='AC'&dkpl$year%in%2008:2010,],FUN=function(x){mean(x,na.rm=T)/2},keep.names=TRUE)
dkps.l$mcfa<-(dkps.l$brt-(dkps.l$lf_br+dkps.l$dBdb))
dkps.l$mcfa.adj<-(dkps.l$brt.adj-(dkps.l$lf_br+dkps.l$dBdb))
dkps.l$mcfr<-1-((dkps.l$lf_br+dkps.l$dBdb)/dkps.l$brt)
dkps.l$mcfr.adj<-1-((dkps.l$lf_br+dkps.l$dBdb)/dkps.l$brt.adj)

dkpsm.l<-summaryBy(NPPbr+brt+lf_br+dBdb+mcfa+mcfr+NPPbr0~site,dkps.l,FUN=me)
dkpsm.adj.l<-summaryBy(NPPbr.adj+brt.adj+lf_br+dBdb+mcfa.adj+mcfr.adj~site,dkps.l,FUN=me)

quartz(w=3.42,h=1.8)
par(mfrow=c(1,1))
par(lwd=.3)
par(mai=c(.4,.4,.2,.4))
plot(NA,xlim=c(0,11.5),ylim=c(-20,350),xaxt='n',xlab='',ylab='', yaxt='n')
mtext('Branch carbon dynamics',line=-.1,adj=0.005,font=1,cex=7/12)
mtext('Duke FACE, 2008-2010',line=-.1,adj=1,font=1,cex=7/12)
mtext('2006-2007',line=.3,adj=1,font=1,cex=7/12)
lines(c(-5,15),c(0,0),lwd=.5)

df.adj.h<-dkpsm.adj.h
df.adj.l<-dkpsm.adj.l
# lines(c(1,4.25),rep(dkpsm.h[1],2),lty=3,lwd=1)
# lines(c(2.5,4.25),rep(dkpsm.h[3],2),lty=3,lwd=1)
# text(5.25,(dkpsm.h[1]+dkpsm.h[3])/2,'Increment',adj=.5,cex=0.4)
# lines(c(1,4.25),rep(dkpsm.l[1],2),lty=3,lwd=1)
# lines(c(2.5,4.25),rep(dkpsm.l[3],2),lty=3,lwd=1)
# text(5.25,(dkpsm.l[1]+dkpsm.l[3])/2,'Increment',adj=.5,cex=0.4)

for(i in seq(1,9,2)){
	polygon(c(-.5,-.5,0,0)+i,c(0,df.adj.h[1,i],df.adj.h[1,i],0),col= 'orange1',border=1,lwd=.75)
	lines(rep((i-0.25),2), df.adj.h[1,i]+c(-1,1)* df.adj.h[1,i+1])
	polygon(c(0,0,.5,.5)+i,c(0,df.adj.l[1,i],df.adj.l[1,i],0),col= 'darkorange2',border=1,lwd=.75)
	lines(rep((i+0.25),2), df.adj.l[1,i]+c(-1,1)* df.adj.l[1,i+1])}

for(i in seq(1,9,2)){	polygon(c(-.5,-.5,0,0)+i,c(0,dkpsm.h[1,i],dkpsm.h[1,i],0),col='green2',border=1,lwd=.75)
	lines(rep((i-0.25),2),dkpsm.h[1,i]+c(-1,1)*dkpsm.h[1,i+1])
polygon(c(0,0,.5,.5)+i,c(0,dkpsm.l[1,i],dkpsm.l[1,i],0),col='green4',border=1,lwd=.75)
	lines(rep((i+0.25),2),dkpsm.l[1,i]+c(-1,1)*dkpsm.l[1,i+1])}

axis (1,seq(0,000,by=2000),tck=.02,label=c(expression(paste('NPP')),expression(paste('Turnover')),expression(paste('F')),expression(paste('∆W'[dead])),'Mass loss'),at=1:5*2-1,mgp=c(0,-.5,0),cex.axis=7/12,lwd=0.3,line=.35)

axis (2,seq(-600,800,by=150),tck=.02,label=TRUE,mgp=c(0,-.1,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Branch compartment')),1,line=.6,font=1,cex=7/12)
mtext(expression(paste('Carbon flux (g C m'^-2,' y'^-1,')')),2,line=.6,font=1,cex=7/12)

adj.tb.h<-dkpsm.adj.h$mcfa.adj.m-dkpsm.h$mcfa.m
ht<-315
polygon(c(-.2,-.2,.2,.2)+5.5,c(ht,ht+adj.tb.h,ht+adj.tb.h,ht),col= 'orange1',border=1,lwd=.5)
polygon(c(-.2,-.2,.2,.2)+5.5,c(ht,ht+adj.tb.h/2,ht+adj.tb.h/2,ht),col= 'darkorange2',border=1,lwd=.5)
text(6,ht,expression(paste('Within-canopy turnover')),adj=c(0,0),cex=7/12)

text(6.6,dkpsm.h$mcfa.m+105,expression(paste('Mass loss / Turnover')),adj=-.2,cex=7/12,font=2)
text(7.2,dkpsm.h$mcfa.m,'50 ± 17%',adj=-1,cex=7/12)
text(7.2,dkpsm.adj.h$mcfa.adj.m,'58 ± 14%',adj=-1,cex=7/12)
text(11.75,dkpsm.l$mcfa.m-5,'-15 ± 20%',adj=1,cex=7/12)
text(11.5,dkpsm.adj.l$mcfa.adj.m+5,'3 ± 2%',adj=1,cex=7/12)





# # ######################################
mod<-read.table('BT_Model.txt',sep='\t',head=TRUE)
mod$Bst<-mod$Vst*400/10000
for (i in 1:2) for (j in 0:1){
mod$NPPst.cum[mod$sp==i&mod$trt==j]<-cumsum(mod$NPPst[mod$sp==i&mod$trt==j])}
mod$rNPPst[mod$trt==0]<-mod$NPPst[mod$trt==0]/mod$NPPst[mod$trt==1]
mod$rNPPst.cum[mod$trt==0]<-mod$NPPst.cum[mod$trt==0]/mod$NPPst.cum[mod$trt==1]
mod$rBst[mod$trt==0]<-mod$Bst[mod$trt==0]/mod$Bst[mod$trt==1]
mod$lty<-ifelse(mod$trt==1,1,4)

range(mod$rNPPst[mod$sp==1&mod$trt==0])
range(mod$rNPPst[mod$sp==2&mod$trt==0])

range(mod$rBst[mod$sp==1&mod$trt==0])
range(mod$rBst[mod$sp==2&mod$trt==0])

ed2_Bst1<-read.table('/Users/hyli0001/Documents/wd/7_Branch_turnover/Model_simulations/ed2_stem_biomass_15Sept.csv',sep=',',head=TRUE)
ed2_Bst<-ed2_Bst1
ed2_Bst$rBst.dk<-ed2_Bst$duke.no.bt/ed2_Bst$duke.bt
ed2_Bst$rBst.fin<-ed2_Bst$finland.no.bt/ed2_Bst$finland.bt
range(ed2_Bst$rBst.dk)
range(ed2_Bst$rBst.fin)

range(mod$Bst[mod$sp==2&mod$trt==1])

tail(ed2_Bst)
tail(mod)

### FIGURES


quartz(w=3.42,h=3.6)
par(mfrow=c(2,2))
par(lwd=.3)
par(mai=c(.05,.35,.55,.1))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,100),ylim=c(0,34.5))
df<-mod[mod$sp==2,]
for (i in 0:1){
	col<-ifelse(i==0,'forestgreen','red3')
	points(Bst~age, df[df$trt==i,],col= col,cex=0.4,lwd=1.5,type='l',lty=1)
	}
mtext('a',line=-.8,adj=0.03,font=2,cex=8/12)
mtext('P. abies',line=-.8,adj=.95,font=3,cex=7/12)
mtext('PREBAS',line=-.1,adj=.5,font=1,cex=7/12)

axis (1,seq(-100,200,by=25),tck=.02,label=FALSE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-100,200,by=10),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Stem biomass (kg m'^-2,')')),2,line=.7,at=-1,font=1,cex=7/12)
legend('bottomright',x.intersp=.5,c(expression(paste('– branch turnover')),expression(paste('+ branch turnover'))),col=c('forestgreen','red3'),lwd=1.5,cex=7/12,box.lwd=0,text.font=3)

par(mai=c(.7,.45,.7,.65),new=TRUE)
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,105),ylim=c(1,1.24))
points(rBst~age,df,col='darkblue',cex=0.3,lwd=1.2,type='l')
axis (1,seq(-100,200,by=50),tck=.02,label=TRUE,mgp=c(0,-.4,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-1,2,by=0.1),tck=.02,label=TRUE,mgp=c(0,-.1,0),cex.axis=7/12,lwd=0.3)

par(mai=c(.05,.05,.55,.4))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,100),ylim=c(0,34.5))	
points(finland.no.bt~year,ed2_Bst,col='forestgreen',cex=0.4,lwd=1.5,type='l',lty=1)
points(finland.bt~year,ed2_Bst,col='red3',cex=0.4,lwd=1.5,type='l',lty=1)

mtext('b',line=-.8,adj=0.03,font=2,cex=8/12)
mtext('P. abies',line=-.8,adj=.95,font=3,cex=7/12)
# mtext('60.9N°, 22.3E',line=-.8,adj=0.95,font=1,cex=5.5/12)
mtext('ED2',line=-.1,adj=.5,font=1,cex=7/12)

axis (1,seq(-100,200,by=25),tck=.02,label=FALSE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-100,200,by=10),tck=.02,label=FALSE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)

par(mai=c(.7,.15,.7,.95),new=TRUE)
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,100),ylim=c(1,1.24))
points(rBst.fin~year,ed2_Bst,col='darkblue',cex=0.3,lwd=1.2,type='l')
axis (1,seq(-100,200,by=50),tck=.02,label=TRUE,mgp=c(0,-.4,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-1,2,by=0.1),tck=.02,label=TRUE,mgp=c(0,-.1,0),cex.axis=7/12,lwd=0.3)

df<-mod[mod$sp==1,]
par(mai=c(.55,.35,.05,.1))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,100),ylim=c(0,72))
for (i in 0:1){
	col<-ifelse(i==0,'forestgreen','red3')
	points(Bst~age, df[df$trt==i,],col=col,cex=0.4,lwd=1.5,type='l',lty=1)
	}
axis (1,seq(-100,200,by=25),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-100,200,by=20),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Age (year)')),1,line=.5,font=1,cex=7/12,at=110)
mtext('c',line=-.8,adj=0.03,font=2,cex=8/12)
mtext('P. sylvestris',line=-.8,adj=.95,font=3,cex=7/12)
# mtext('60.9N°, 22.3E',line=-.8,adj=0.98,font=1,cex=5.5/12)
# mtext('PREBAS',line=-.1,adj=1,font=1,cex=6/12)


par(mai=c(1.1,.45,.3,.65),new=TRUE)
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,100),ylim=c(1,1.24))
points(rBst~age,df,col='darkblue',cex=0.3,lwd=1.2,type='l')
axis (1,seq(-100,200,by=50),tck=.02,label=TRUE,mgp=c(0,-.4,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-1,2,by=0.1),tck=.02,label=TRUE,mgp=c(0,-.1,0),cex.axis=7/12,lwd=0.3)
# mtext(expression(paste('Ratio')),4,line=-.3,font=1,cex=5/12)


par(mai=c(.55,.05,.05,.4))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,100),ylim=c(0,72))
points(duke.no.bt~year,ed2_Bst,col='forestgreen',cex=0.4,lwd=1.5,type='l',lty=1)
points(duke.bt~year,ed2_Bst,col='red3',cex=0.4,lwd=1.5,type='l',lty=1)
mtext('d',line=-.8,adj=0.03,font=2,cex=8/12)
mtext('P. taeda',line=-.8,adj=.95,font=3,cex=7/12)
# mtext('36.0N°, 79.0W',line=-.8,adj=0.85,font=1,cex=5.5/12)
axis (1,seq(-100,200,by=25),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-100,200,by=20),tck=.02,label=FALSE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
# mtext('ED2',line=-.1,adj=1,font=1,cex=6/12)


par(mai=c(.75,.65,.65,.45),new=TRUE)
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,100),ylim=c(1,1.24))
points(rBst.dk~year,ed2_Bst,col='darkblue',cex=0.3,lwd=1.2,type='l')
axis (1,seq(-100,200,by=50),tck=.02,label=TRUE,mgp=c(0,-.4,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-1,2,by=0.1),tck=.02,label=TRUE,mgp=c(0,-.1,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('ratio')),2,line=.35,font=1,cex=6/12)
mtext(expression(paste('age')),1,line=.05,font=1,cex=6/12)















# ### Supplementary FIGURES
	# df.of<-smg[smg$site%in%c("Bräc","Ebbe","Grän","Gävl","Möln"),]

# # quartz(w=3.5,h=3.2)
# par(mfrow=c(2,2))
# par(lwd=.3)
# par(mai=c(.2,.5,.4,.1))
# plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0.2,1),ylim=c(0,.1))
# points(brtr~h.d, df.of,col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
# mtext('a',line=-.1,adj=0.01,font=2,cex=8/12)
# summary(lm(h.d.pre~brtr,smg[smg$sp==1,]))
# # mtext(expression(paste('r'^2,'=0.90',sep='')),line=-.8,adj=.95,font=1,cex=5/12)
# # mtext(expression(paste('p<.001',sep='')),line=-1.3,adj=.95,font=1,cex=5/12)

# axis (1,seq(-2,9,by=2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
# axis (2,seq(-1,2,by=.5),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
# mtext(expression(paste('Branch turnover rate')),2,line=1.2,font=1,cex=6/12)
# mtext(expression(paste('(kg T'[BR],' kg'^-1,' W'[BR],' yr'^-1,')')),2,line=.5,font=1,cex=6/12)
# mtext(expression(paste('Height increment (∆H; m yr'^-1,')')),1,line=.5,font=1,cex=6/12)

	# site<-unique(df.of$site[df.of$site!='hb'])
	# for (j in site){
		# df1<-df.of[df.of$site==j,]
		# ttrt<-unique(df1$ttrt)
		# for (k in ttrt){
			# df2<-df1[df1$ttrt==k,]
			# df3<-summaryBy(brtr+h.d~bg+pch+col,df2,FUN=me)
			# lines(rep(df3$h.d.m,2),c(-1,1)*df3$brtr.se+df3$brtr.m)
			# lines(c(-1,1)*df3$h.d.se+df3$h.d.m,rep(df3$brtr.m,2))
			# points(brtr.m~h.d.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=0.2)
			# }}
		# points(brtr~h.d, df.of[df.of $site=='hb',],bg=bg,pch=pch,col=col,cex=7/12,lwd=0.2)
	# curve(h.d(x,1),add=TRUE,xlim=range(smg$h.d[smg$sp==1],na.rm=TRUE),lwd=0.5)

# lnd1<-rbind(lnd[lnd$sp==1,][1:14,])
# legend('topleft',x.intersp=.5,y.intersp=.85,lnd1$st,pch=21,cex=0.4,col=lnd1$col,pt.bg=lnd1$bg,box.col=0,ncol=5,pt.lwd=.1)
# lnd2<-rbind(lnd4,lnd[lnd$sp==2,][c(1:6,8),],lnd4,lnd[lnd$sp==2,][7,],lnd5)
# legend('topleft',x.intersp=.5,y.intersp=.85,lnd2$st,pch=21,cex=0.4,col=lnd2$col,pt.bg=lnd2$bg,box.col=0,ncol=2,pt.lwd=.75)


