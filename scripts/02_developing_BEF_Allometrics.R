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




exp_decay<- function(x, L, A, k) {
  L + A * exp(-k * x)
}

mich_men <- function(x,L,A,r0) {
  L + A / (1 + x/r0)
}


library(brms)
library(cmdstanr)
fit_bef_joint_bayes <- function(data,
                                y_1,
                                y_2,
                                x,
                                group,
                                chains = 4,
                                iter = 4000,
                                cores = 4,
                                adapt_delta = 0.95,
                                max_treedepth = 12) {
	df <- data[, c(y_1, y_2, x, group)]
	names(df) <- c("y1", "y2", "x", "grp")

	df <- df[complete.cases(df[, c("x", "grp")]), ]
	df$grp <- factor(df$grp)
	med_x <- median(df$x, na.rm = TRUE)
	prior_r0 <- paste0("lognormal(", log(med_x), ", 1)")
	
	bf_y1 <- bf(
		y1 ~ L + A / (1 + x / r0),
		L  ~ 1 + (1 | grp),
		A  ~ 1 + (1 | grp),
		r0 ~ 1 + (1 | grp),
		nl = TRUE)

	bf_y2 <- bf(
		y2 | mi() ~ L + A / (1 + x / r0),
		L  ~ 1 + (1 | grp),
		A  ~ 1 + (1 | grp),
		r0 ~ 1 + (1 | grp),
		nl = TRUE)

fit <- brm(
	bf_y1 + bf_y2 + set_rescor(FALSE),
	data = df,
	family = gaussian(),
	prior = c(
		prior(normal(1.2, 0.5), nlpar = "L",  resp = "y1", lb = 0),
		prior(normal(0.6, 0.6), nlpar = "A",  resp = "y1", lb = 0),
		prior_string(prior_r0, nlpar = "r0", resp = "y1", lb = 0),

		prior(normal(1.4, 0.5), nlpar = "L",  resp = "y2", lb = 0),
		prior(normal(0.7, 0.7), nlpar = "A",  resp = "y2", lb = 0),
		prior_string(prior_r0, nlpar = "r0", resp = "y2", lb = 0)),
	chains = chains,
	iter = iter,
	cores = cores,
	control = list(
		adapt_delta = adapt_delta,
		max_treedepth = max_treedepth),
	backend = "cmdstanr",
	file = "fit_bef_joint",
	file_refit = "on_change",
	save_pars = save_pars(all = TRUE)
	)
	return(fit)
}

fit_joint<-fit_bef_joint_bayes(bp,y_1="befa.st",y_2="beft.st",x="sdi",group = "ft1.forest_type",2,400,2)
plot(fit_joint)
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

non.mod<-nls(befa.st~L+A/(1+sdi/r0),bp[bp$Family=='Pinaceae',],start=c(L=0.3,A=1.2,r0=2.6))

plot(befa.st~sdi,bp[bp$Family=='Pinaceae',],col=0)
for(i in 1:6){
	sp_c<-unique(bp[bp$Family=='Pinaceae','sp_code'])[i]
points(befa.st~sdi,bp[bp$sp_code==sp_c&bp$Family=='Pinaceae',],bg=i,pch=21)}
points(befa.st~sdi,bp[bp$sp_code=='LL'&bp$Family=='Pinaceae',],bg='white',pch=22,col=4)
df00<-bp[bp$Family=='Pinaceae',]
sdi<-as.data.frame(seq(range(df00$sdi)[1],range(df00$sdi)[2],1))
colnames(sdi)<-'sdi'
sdi$befa.a.pre<-predict(non.mod,sdi)
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
