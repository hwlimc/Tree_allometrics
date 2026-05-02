# renv::restore('..')
# renv::status('..')
getwd()
source('01_data_process_plot_havests.R')
bp<-read.table('../processed_data/plot_biomass.txt',sep='\t',head=TRUE)
bp$d2h<-bp$d^2*bp$h
bp$cai<-(bp$Vst-bp$Vst.5)/5
bp$Bw<-(bp$Bst+bp$Bbr+bp$Bcr)
bp$Baw<-(bp$Bst+bp$Bbr)
bp$Ba<-(bp$Bst+bp$Bbr+bp$Bf)
bp$Bt<-(bp$Bst+bp$Bbr+bp$Bcr+bp$Bf)

bp$Bst.5<-bp$Vst.5*bp$bwd
bp$Bbr.5<-bp$Bst.5*(bp$Bbr/bp$Bst)
bp$Bf.5<-bp$Bst.5*(bp$Bf/bp$Bst)
bp$Bcr.5<-bp$Bst.5*(bp$Bcr/bp$Bst)

bp$Baw.5<-(bp$Bst.5+bp$Bbr.5)
bp$Bw.5<-(bp$Baw.5+bp$Bcr.5)
bp$Ba.5<-(bp$Baw.5+bp$Bf.5)
bp$Bt.5<-(bp$Ba.5+bp$Bcr.5)

bp$dBst<-(bp$Bst-bp$Bst.5)/5
bp$dBw<-(bp$Bw-bp$Bw.5)/5
bp$dBaw<-(bp$Baw-bp$Baw.5)/5
bp$dBt<-(bp$Bt-bp$Bt.5)/5
bp$mBf<-(bp$Bf+bp$Bf.5)/2
bp$ge.aw<-bp$dBaw/bp$mBf
bp$ge.w<-bp$dBw/bp$mBf
bp$ge.ab[bp$PTF%in%c('DBF','DNF')]<-(bp$dBaw[bp$PTF%in%c('DBF','DNF')]+bp$mBf[bp$PTF%in%c('DBF','DNF')])/bp$mBf[bp$PTF%in%c('DBF','DNF')]
bp$ge.t[bp$PTF%in%c('DBF','DNF')]<-(bp$dBw[bp$PTF%in%c('DBF','DNF')]+bp$mBf[bp$PTF%in%c('DBF','DNF')])/bp$mBf[bp$PTF%in%c('DBF','DNF')]


par(mfrow=c(2,3))
for (i in unique(bp$sp_code)){
	df1<-bp[bp$sp_code==i,]
	plot(befa.st~ species.richness,df1,main=i,ylim=c(1,5))
	points(beft.st~ species.richness,df1,col=2)

	plot(befa.st~ species.shannon,df1,ylim=c(1,5))
	points(beft.st~ species.shannon,df1,col=2)

	plot(befa.st~ cv,df1,ylim=c(1,5))
	points(beft.st~ cv,df1,col=2)

	plot(befa.st~ wb.scale,df1,ylim=c(1,5))
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,df1,ylim=c(1,5))
	points(beft.st~ sdi,df1,col=2)
	
	plot(befa.st~dist.mmd,df1,ylim=c(1,5))
	points(beft.st~ dist.mmd,df1,col=2)	
	
	}
head(df1)
par(mfrow=c(2,3))
for (i in unique(bp$PFT[!is.na(bp$Vst.5)])){
	df1<-bp[bp$PFT==i&!is.na(bp$Vst.5),]
	plot(ge.aw~PFT1.richness,df1,main=i)
	points(ge.w~PFT1.richness,df1,col=2)

	plot(ge.aw ~ PFT1.shannon,df1)
	points(ge.w ~ PFT1.shannon,df1,col=2)

	plot(ge.aw ~ cv,df1)
	points(ge.w ~ cv,df1,col=2)

	plot(ge.aw ~ wb.scale,df1)
	points(ge.w ~ wb.scale,df1,col=2)

	plot(ge.aw ~sdi,df1)
	points(ge.w ~ sdi,df1,col=2)
	
	plot(ge.aw ~dist.mmd,df1)
	points(ge.w ~ dist.mmd,df1,col=2)		
	}



