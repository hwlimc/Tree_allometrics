library("doBy")
library("segmented")
library("nlme")
library("lme4")
library("reshape")
#library('rworldmap')
#library('rworldxtra')
#library('Evapotranspiration')
library('emmeans')
library('multcomp')
library('plyr')
library('agricolae')
library('multcompView')
library('rmarkdown')
# library('colorout')
# setOutputColors(normal = 4, negnum = 4, zero = 4, number = 4,
                     # date = 1, string = 6, const = 5, false = 5,
                     # true = 2, infinite = 5, index = 0, stderror = 4,
                     # warn = c(1, 0, 1), error = c(1, 7),
                     # verbose = TRUE, zero.limit = NA)

# library('rethinking')
# library("rstan")
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

## cast(data, col.group1+col.group2~row.group,FUN)
rm(list=ls())
'%out%'<-function(a,b){!a%in%b}
se<-function(x){sd(x,na.rm=T)/sqrt(sum(!is.na(x)))}
cv<-function(x){se(x)/mean(x,na.rm=T)}
n<-function(x){sum(!is.na(x))}
md<-function(x){c(m=mean(x,na.rm=T),sd=sd(x,na.rm=T))}
me<-function(x){c(m=mean(x,na.rm=T),se=se(x))}
men<-function(x){c(m=mean(x,na.rm=T),se=se(x),n=n(x))}
mn<-function(x){c(m=mean(x,na.rm=T),n=n(x))}
sn<-function(x){c(s=sum(x,na.rm=T),n=n(x))}
rmse<-function(x){sqrt(mean(x^2))}
mse.model<-function(observation,prediction){
	sum((prediction-observation)^2)/sum(observation^2)}
gini<-function(x){ineq(x,type="Gini")}
maxn<-function(x,n){mean(x[order(-x)][1:n])}
max.prop<-function(x,proportion){
	n.prop<-ceiling(sum(!is.na(x))*proportion)
	x<-x[order(-x)]
	mean(x[1:n.prop],na.rm=TRUE)}
min.prop<-function(x,proportion){
	n.prop<-ceiling(sum(!is.na(x))*proportion)
	x<-x[order(x)]
	mean(x[1:n.prop],na.rm=TRUE)}

rep.df<-function(df,rep.time){
	df[rep(1:nrow(df),rep.time),]}
	
pred.fit<-function(model,df,x,col){
	mn<-min(df[,x],na.rm=TRUE)
	mx<-max(df[,x],na.rm=TRUE)
	pred<-seq(mn,mx,length.out=300)
	response<-as.data.frame(pred)
	colnames(response)<-x
	m<-predict(model,response)
	lines(pred,m,col=col)
	}

pred.err<-function(model,df,x,col){
	mn<-min(df[,x],na.rm=TRUE)
	mx<-max(df[,x],na.rm=TRUE)
	pred<-seq(mn,mx,length.out=300)
	response<-as.data.frame(pred)
	colnames(response)<-x
	m<-predict(model,response)
	err<-summary(model)$sig
	nse<-m-(err*2)
	pse<-m+(err*2)
	polygon(c(pred,rev(pred)),c(nse,rev(pse)),col=col,border=FALSE)
	}


mcolname<-function(df,nr.group,remove,add){c(colnames(df)[1:(nr.group)],paste(substr(colnames(df)[(nr.group+1):length(df)],1,nchar(colnames(df)[(nr.group+1):length(df)])-remove),add,sep=""))}


# coeff.list <- dlply(NPP, .(year,teat), function(DF) lm(NPP~sd, data=DF))
# coef<-ldply(coeff.list, extractfun)
extractfun <- function(m) {
       cf <- coef(m)
       tinfo <- summary(m)$coefficients[2, c(2, 4)]
       r2 <- summary(m)$r.squared
       data.frame(intercept = cf[1], slope = cf[2],
                  slope.se = tinfo[1], n = length(resid(m)),
                  pval = tinfo[2], Rsq = r2)}

#between variables (#assuming unequal variance)

aovBy <- function(formula, group, data, ...){
    formulaFunBy(formula, group, data, FUN=summary, class="aovBy", ...)
    }

t.testBy1 <- function(formula, group, data, ...){
    formulaFunBy(formula, group, data, FUN=t.test, class="t.testBy1", ...)
    }
    
#between columns
t.testBy2 <- function(formula, group, data, ...){
    xyFunBy(formula, group, data, FUN=t.test, class="t.testBy1", ...)
  }

#t-test with mean and standard error
#Equality of variance
two.s.f<-function(s1,s2,n1,n2){
	sd1<-s1*sqrt(n1)
	sd2<-s2*sqrt(n2)
	2*ifelse(sd1>sd2,(1-pf((sd1^2*(n1))/(sd2^2*(n1)),n1-1,n2-1)),(1-pf((sd2^2*(n2))/(sd1^2*(n1)),n2-1,n2-1)))
	}
	
#Independent (unpaired) one-sample test (equal or unequal sample size, and equal variance)
one.s.t<-function(x,a,s,n){
	(1-pt(abs((x - a)/s), n-1))*2
	}
	


#Independent (unpaired) two-sample test (equal or unequal sample size, and equal variance)
two.s.t<-function(x1,x2,s1,s2,n1,n2){
	(1-pt(abs((x1-x2)/(sqrt(((n1-1)*(s1*sqrt(n1))^2+(n2-1)*(s2*sqrt(n2))^2)/(n1+n2-2))*sqrt(1/n1+1/n2))), n1+n2-2))*2
	}
#Independent (unpaired) two-sample test (equal or unequal sample size, and unequal variance), Welch's t-test
Wel.df<-function(x1,x2,s1,s2,n1,n2){(s1^2+s2^2)^2/(s1^4/(n1-1)+s2^4/(n2-1))}
two.s.t.uv<-function(x1,x2,s1,s2,n1,n2){
	2*(1-pt(abs((x1-x2)/(sqrt(s1^2+s2^2))),Wel.df(x1,x2,s1,s2,n1,n2)))
	}
#Independent (unpaired) two-sample test (equal or unequal sample size, and combination of unequal and equal variances)	
two.tftest<-function(x1,x2,s1,s2,n1,n2){
	fnc<-function(x1,x2,s1,s2,n1,n2){
	ftest<-two.s.f(s1,s2,n1,n2)
	ifelse(ftest>0.05,two.s.t(x1,x2,s1,s2,n1,n2),two.s.t.uv(x1,x2,s1,s2,n1,n2))}
	mapply(fnc,x1,x2,s1,s2,n1,n2)}


#Power calculation #NEED TO BE FIXED
pwr.t<-function(x1,x2,s1,s2,n1,n2){
a <- x2-x1
s <- sqrt(s1*s1+s2*s2)
ftest<-two.s.f(s1,s2,n1,n2)
df<-ifelse(ftest>0.05,n1+n2-2,Wel.df(x1,x2,s1,s2,n1,n2))
error <- qt(0.975,df=df)*s
left <- a-error
right <- a+error
tleft <- left/s
tright <- right/s
1-(pt(tright,df=df)-pt(tleft,df=df))}


sd.mean<-function(sd){
	sqrt(mean(sd^2,na.rm=TRUE))}
sd.m<-function(m1,sd1,m2){
	sqrt((m2/m1)^2)*sd1}
sd.pooled<-function(sd1,n1,sd2,n2){
	sqrt(((n2-1)*sd2^2+((n1-1)*sd1^2))/(n2+n1-2))}
sd.mtom<-function(m1,sd1,m2,sd2){
	sqrt((sd1/m1)^2+(sd2/m2)^2)}
sd.mtom1<-function(m,sd){
	sqrt(sum((sd/m)^2))}	
sd.sum<-function(sd){
	sqrt(sum(sd^2))}	# x = SD or SE
sd.msubm<-function(sd1,sd2){
	sqrt(sd1^2+sd2^2)}
sd.mtob<-function(a,m,b,sd){
	abs((a*m^b*b*sd)/m)}
log.sd<-function(a,m,sd,lower){abs(a*(sd/(m*log(lower))))}

effect.size<-function(m1,m2,sd1,sd2,n1,n2){
	(m2-m1)/sd.pooled(sd1,n1,sd2,n2)}

#Stand Density Index
SDI<-function(N.ha,D.cm,coef){
	((N.ha/2.4711)*(D.cm/25.4)^coef)}
qdr.m<-function(D){
	sqrt(mean(D^2,na.rm=TRUE))}
BA<-function(D){
	D^2/4*pi}
revD<-function(BA){
	sqrt((BA*4)/pi)}


#VPD kPa
vpd<-function(t,rh){
	exp(t*17.269/(t+237.3)*.61078)*(1-(rh/100))}

vpdtxn<-function(Tx,Tn){((6.1078 * exp(17.269 * Tx / (237.3 + Tx)))-(6.1078 * exp(17.269 * Tn / (237.3 + Tn))))/2}	#mBar

vpd1<-function(T,rh){
A<--1.88*10^4
B<--13.1
C<--1.5*10^-2
D<-8*10^-7
E<--1.69*10^-11
F<-6.456
exp(A/(T+273.15)+B+C*(T+273.15)+D*(T+273.15)^2+E*(T+273.15)^3+F*log(T+273.15))-exp(A/(T+273.15)+B+C*(T+273.15)+D*(T+273.15)^2+E*(T+273.15)^3+F*log(T+273.15))*rh/100}

get.vpd <- function(temp,rh){
  ## calculate vapor pressure deficit
  vpd <- ((100 - rh) / 100) * 6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 / (273 + temp)))
  return(vpd)
}

## Stem Volume functions for pines
NäshlundV.pine<-function(d.m,h.m){
	#return in m3
	#for diameter < 5 cm
	d<-d.m*100
	h<-h.m
	ifelse(is.na(d)| is.na(h),NA,ifelse(d>=5,NaN,(0.22+0.1066*(d^2)+0.02085*(d^2*h)+0.008427*d*h^2)/1000))}

BrandelV.pine<-function(d.m,h.m){
	#return in m3
	#for diameter ≥ 5cm
	d<-d.m*100
	h<-h.m
	ifelse(is.na(d)| is.na(h),NA,ifelse(d<5 | h<=1.3,NaN,(10^(-1.38903+1.84493*log(d,10)+0.06563*log(d+20,10)+2.02122*log(h,10)-1.01095*log(h-1.3,10)))/1000))}

BrandelV1.pine<-function(d.m,h.m){
	#return in m3
	#for diameter ≥ 5cm
	d<-d.m*100
	h<-h.m
	ifelse(is.na(d)| is.na(h),NA,ifelse(d<5 | h<=1.3,NaN,(10^-1.38903*d^1.84493*(d+20)^0.06563*h^2.02122*(h-1.3)^-1.01095)/1000))}

Vst.estimation.pine<-function(d.m,h.m){
	#return in m3
	ifelse(d.m<0.05,	NäshlundV.pine(d.m,h.m),BrandelV.pine(d.m,h.m))}


#write.table(NPP,"NPP.txt",sep="\t",quote=FALSE,row.names=FALSE)
#plot(getMap(resolution="high"),xlim=c(15,20),ylim=c(55,70))
	
# lwd 1 = 0.75 pt
# cex 1 = 12 pt
