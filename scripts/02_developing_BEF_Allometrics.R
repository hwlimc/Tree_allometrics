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


exp_decay<- function(RSD, L, A, k) {
  L + A * exp(-k * RSD)
}

fit_bef <- function(df) {
  nls(
    BEF ~ L + A * exp(-k * RSD),
    data = df,
    start = list(
      L = min(df$BEF, na.rm = TRUE),
      A = max(df$BEF, na.rm = TRUE) - min(df$BEF, na.rm = TRUE),
      k = 1
    )
  )
}

init_k <- coef(lm(log(BEF) ~ RSD, data = df))[2] * -1

fit <- nls(
  BEF ~ L + A * exp(-k * RSD),
  data = df,
  start = list(
    L = min(df$BEF),
    A = diff(range(df$BEF)),
    k = abs(init_k)
  )
)

mich_men <- function(rsd, bef0, amp, r0) {
  bef0 + amp / (1 + rsd/r0)
}

fit_befa.st <- function(df) {
	nls(befa.st~bef0 + amp / (1 + sdi/r0),
		data = df,
	    start = list(bef0 = min(df$BEF, na.rm = TRUE),
						amp = max(df$BEF, na.rm = TRUE) - min(df$BEF, na.rm = TRUE),
						r0 = median(df$RSD,na.rm=TRUE)),
		algorithm = "port",
		lower = c(bef0 = 0, amp = 0, r0 = 1e-10)
		)
}

df1<-bp[bp$ft1.forest_type =='mono_N',]

mich_ment.mod<- function(df,y,x) {
	df.mod<-data.frame(xx=df[[x]],yy=df[[y]])
	df.mod<-df.mod[complete.cases(df.mod),]
	nls(yy~L + A / (1 + xx/r0),
		data = df.mod,
	    start = list(L = min(df.mod$yy, na.rm = TRUE),
	    			A = max(df.mod$yy, na.rm = TRUE) - min(df.mod$yy, na.rm = TRUE),
					r0 = median(df.mod$xx,na.rm=TRUE)),
		algorithm = "port",
		lower = c(bef0 = 0, amp = 0, r0 = 1e-10)
		)
}
fit<-mich_ment.mod(df1,'befa.st','sdi')
	    
newx<-seq(min(df1$sdi, na.rm=TRUE),max(df1$sdi, na.rm=TRUE),length.out=100)
pred<-predict(fit,newdata = data.frame(xx = newx))

	plot(befa.st~sdi,ylim=c(1,lmx),df1)
	points(beft.st~ sdi,df1,col=2)
	lines(newx, pred, col=4, lwd=2)

par(mfrow=c(2,2))
for (i in unique(bp$ft1.forest_type)){
	df1<-bp[bp$ft1.forest_type ==i,]
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
