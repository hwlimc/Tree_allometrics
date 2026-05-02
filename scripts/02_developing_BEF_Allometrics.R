# renv::restore('..')
# renv::status('..')
getwd()
# source('01_data_process_plot_havests.R')
# system('ls -alt ../processed_data')
bp<-read.table('../processed_data/plot_biomass.txt',sep='\t',head=TRUE)
bp$d2h<-bp$d^2*bp$h
bp$cai<-(bp$Vst-bp$Vst.5)/5
bp$Bw<-(bp$Bst+bp$Bbr+bp$Bcr)
bp$Baw<-(bp$Bst+bp$Bbr)
bp$Ba<-(bp$Bst+bp$Bbr+bp$Bf)
bp$Bt<-(bp$Bst+bp$Bbr+bp$Bcr+bp$Bf)
sum(!is.na(unique(bp$sp_code)))
unique(bp$sp_code)
unique(bp$ft1.forest_type)

par(mfrow=c(2,2))
for (i in unique(bp$ft1.forest_type)){
	df1<-bp[bp$ft1.forest_type ==i,]
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

par(mfrow=c(2,2))
for (i in unique(bp$Family)){
	df1<-bp[bp$Family==i,]
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
	

par(mfrow=c(2,2))
for (i in unique(bp$Genus)){
	df1<-bp[bp$Genus==i,]
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


par(mfrow=c(2,2))
for (i in unique(bp$sp_code)[1:12]){
	df1<-bp[bp$sp_code==i,]
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
