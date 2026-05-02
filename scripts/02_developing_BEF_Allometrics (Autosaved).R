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

unique(bp$PFT)
bp$FoT<-ifelse(bp$dom.sp.prop>0.75,ifelse(substr(bp$PFT,1,1)=='E','D')

plot(bp$dom.sp.prop)

par(mfrow=c(2,2))
for (i in unique(bp$sp_code)){
	df1<-bp[bp$sp_code==i,]

	plot(befa.st~ cv,df1,ylim=c(1,5))
	points(beft.st~ cv,df1,col=2)

	plot(befa.st~ wb.scale,df1,ylim=c(1,5))
	points(beft.st~ wb.scale,df1,col=2)

	plot(befa.st~sdi,df1,ylim=c(1,5))
	points(beft.st~ sdi,df1,col=2)
	
	plot(befa.st~dist.mmd,df1,ylim=c(1,5))
	points(beft.st~ dist.mmd,df1,col=2)	
	
	}

par(mfrow=c(2,3))
for (i in unique(bp$PFT)){
	df1<-bp[bp$PFT==i,]

	plot(ge.aw ~ cv,df1)
	points(ge.w ~ cv,df1,col=2)

	plot(ge.aw ~ wb.scale,df1)
	points(ge.w ~ wb.scale,df1,col=2)

	plot(ge.aw ~sdi,df1)
	points(ge.w ~ sdi,df1,col=2)
	
	plot(ge.aw ~dist.mmd,df1)
	points(ge.w ~ dist.mmd,df1,col=2)		
	}



