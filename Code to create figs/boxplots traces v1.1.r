basefile	<-	'~/Dropbox/Shared Vito/cFP paper/Figures/'
#Code to automatically add PC compatibility
basefile2	<-	gsub('\\\\','/',path.expand(basefile))
xx			<-	try(setwd(basefile2))
if(class(xx)=="try-error"){basefile2 <- paste('C:/Users/Peter',substring(basefile,2),sep="")}
read.dir 	<- paste(basefile2,"Data used for figures (pulled by R)",sep="")
write.dir	<- paste(basefile2,"Figure parts",sep="")
setwd(read.dir)

##################################################################
# Pull all data
##################################################################

source('cFP norm (72h).R')
source('cFP combo estimate-slopes.r')
library(scales) #for alpha blending
library(gplots) #for error bars: plotCI

##################################################################
#Subset the data

E.TR500		<- subset(cfp, Erl == 3 & TRM == 500)
D.TR500		<- subset(cfp, Erl == 0 & TRM == 500)
E.CHX500	<- subset(cfp, Erl == 3 & CHX == 500)
D.CHX500	<- subset(cfp, Erl == 0 & CHX == 500)
D.CHX50		<- subset(cfp, Erl == 0 & CHX == 50)
D.CHX5		<- subset(cfp, Erl == 0 & CHX == 5)
D.17A		<- subset(cfp, Erl == 0 & X17A == 0.1)
D.ctrl		<- subset(cfp, Erl == 0 & CHX == 0 & FSK == 0 & TRM == 0 & X17A == 0)
IDCHX500	<- unique(D.CHX500$id)

DCHX.slope		<- subset(slopes,Erl==0 & CHX==500)
DCHX.slope50	<- subset(slopes,Erl==0 & CHX==50)
DCHX.slope5		<- subset(slopes,Erl==0 & CHX==5)
D.slope			<- subset(slopes, Erl == 0 & CHX == 0 & FSK == 0 & TRM == 0 & X17A == 0)

#####################################################################
#Plot average dynamics in 500 nM CHX versus control: Fig. 3a

ctrl.mean	<-	aggregate(D.ctrl$nl2,by=list(D.ctrl$Time.day),FUN=mean)
colnames(ctrl.mean)	<-	c("Time.day","nl2")
CHX.mean	<-	cbind(
					aggregate(D.CHX500$nl2,by=list(D.CHX500$Time.day),FUN=mean),
					aggregate(D.CHX500$nl2, by=list(D.CHX500$Time.day),FUN=sd)$x
					)
colnames(CHX.mean)	<-	c("Time.day","nl2","sd")

dev.new(width=3,height=4)
fn2	<-	paste(write.dir,"/Fig 3 population average dynamics.pdf",sep="")
#pdf(fn2,width=3, height=4)
par(font.lab=2)
plot(CHX.mean$Time.day,CHX.mean$nl2,type='o',ylim=c(0,2),xlab="Time in drug (d)",ylab="Population doublings",lwd=2,pch=26)
lines(ctrl.mean$Time.day,ctrl.mean$nl2,lty=2,lwd=2)

plotCI(
	x	=	CHX.mean$Time.day,
	y	=	CHX.mean$nl2,
	uiw	=	CHX.mean$sd, 
	pch	=	20,
	gap	=	0,
	add	=	TRUE,
	col	=	alpha('black',0.5)
)
legend("topright",c("Ctrl","Chx"),lty=c(1,2),lwd=c(2,2),bty="n")

#rect(3,-0.2,10.3,2.2,col=alpha('lightblue',0.5),border=NA)
#text(x=5,y=0.1,"Linear")
legend("topright",c("Ctrl","Chx"),lty=c(1,2),lwd=c(2,2),bty="n")
#dev.off()
##################################################################
#Plot traces of individual colonies over colony size boxplots: Fig. 3b

lo	<-	"B04_J_E0_CHX500_TRM0_X17A0_Plate1"
hi	<-	"B03_A_E0_CHX500_TRM0_X17A0_Plate1"
down<-	"B05_H_E0_CHX500_TRM0_X17A0_Plate1"
dev.new(width=3,height=4)	
fn5	<-	paste(write.dir,"/Fig 3 boxplot traces.pdf",sep="")
#pdf(fn5,width=3, height=4)
par(cex.lab=1.2)	
boxplot(nl2 ~ Time.day, data=D.CHX500, notch=TRUE, xlab="Time in drug (d)", ylab="Colony doublings", 
		border='grey50')
lines(1:11, subset(cfp,id==lo)$nl2[1:11], col='blue', lwd=2)
lines(1:11, subset(cfp,id==hi)$nl2[1:11], col='red', lwd=2) 
lines(1:11, subset(cfp,id==down)$nl2[1:11], col='green', lwd=2)
lines(1:4, mean.ctrl,col='black',lty=2,lwd=2)

##################################################################
#72h normalization of plot traces of individual colonies over colony size boxplots: Fig. 3d

dev.new(width=3, height=4)
fn.boxplot	<-	paste(write.dir,"/Fig 3 Colony dynamics boxplots.pdf",sep="")
#pdf(fn.boxplot,width=3, height=4)
par(font.lab=2)
boxplot(nl2.72 ~ Time.day, data=D.CHX500[D.CHX500$Time.day>2,], notch=TRUE, xlab="Time in drug(d)", ylab="Colony Doublings", 
		border='grey20',ylim=c(-1,2))
n	<- length(unique(D.CHX500$id))
text(y=2,x=6,paste("n = ",n,sep=""))

lines(1:8, subset(cfp,id==lo)$nl2.72[4:11], col='blue', lwd=2)
lines(1:8, subset(cfp,id==hi)$nl2.72[4:11], col='red', lwd=2)
lines(1:8, subset(cfp,id==down)$nl2.72[4:11], col='green', lwd=2)
mean.ctrl	<-	rep(0,4)
for (i in 1:4) {mean.ctrl[i]<-mean(subset(D.ctrl,Time.day==i-1)$nl2.72)}
mean.ctrl	<-	mean.ctrl-mean.ctrl[1]
lines(1:4, mean.ctrl,col='black',lty=2,lwd=2)
#dev.off()

##################################################################
#Code to display estimated data from linear model fits: Fig. 3e

r	<-	2 #registration error of boxplots		
dev.new(width=3,height=4)
fn4	<-	paste(write.dir,"/Fig 3 linear model colony dynamics.pdf",sep="")
#pdf(fn4,width=3, height=4)
par(font.lab=2)
plot(	4:11-r, 4:11-r, type='n', 
		xlab="Time in drug (d)", ylab="Colony Doublings",
		xlim=c(2.5,10.5)-r,ylim=c(-1,2),xaxt='n')
for (i in 1:nrow(DCHX.slope))
{
	segments(3-r,0,10-r,7*DCHX.slope$Slope[i],col='grey75')
}
lm.lo	<-	coef(lm(subset(cfp,id==lo&Time.day>2)$nl2.72~subset(cfp,id==lo&Time.day>2)$Time.day))[2]
lm.hi	<-	coef(lm(subset(cfp,id==hi&Time.day>2)$nl2.72~subset(cfp,id==hi&Time.day>2)$Time.day))[2]
sn		<-	s[4:11]-s[4]
lm.mean	<-	coef(lm.mean<-lm(sn~c(3:10)))[2]

est.col	<-	matrix(NA,nrow=0,ncol=2)
colnames(est.col)	<-	c("Col.size","Time.day")
for (i in c(1:8)){
	x1	<-	DCHX.slope$Slope*(i-1)		
	#x2	<-	rep(i+2,n)
	x2	<-	rep(i-1,n)
	est.col		<-	rbind(est.col,cbind(x1,x2))
}

boxplot(Col.size ~ Time.day,data=est.col,add=TRUE,xaxt='n',notch=TRUE,border="grey20")
axis(1,at=1:8,labels=c(3:10))
segments(3-r,0,10-r,7*lm.hi,col='red',lwd=2)
segments(3-r,0,10-r,7*lm.lo,col='blue',lwd=2)
segments(3-r,0,10-r,7*lm.mean,col='black',lwd=2)
segments(3-r,0,6-r, 3*mean(D.slope$Slope),col='black',lty=2,lwd=2)

##################################################################
#Code to plot the linear model fit r-squared values for CHX-treated colonies: Supplementary Figure

dev.new(width=3,height=4)
fn8	<-	paste(write.dir,"/Fig 4 R2 linear model coefficients.pdf",sep="")
#pdf(fn8,width=3, height=4)
par(font.lab=2)
R2	<-	append(D.slope$R2,DCHX.slope$R2)
#R2	<-	append(	append(	append(D.slope$R2, DCHX.slope$R2),DCHX.slope50$R2),DCHX.slope5$R2)
hist(R2,xlab="R2",main="",breaks=20)
text(x=0,y=250,paste("n = ",length(R2),sep=""))
#dev.off()

##################################################################
#Code to generate a correlation plot of rate and 10-day change in colony size (outcome): Supplementary Figure

ten.day	<-	subset(D.CHX500,Time.day==10)
pred	<-	cbind(ten.day[,c('id','nl2')],DCHX.slope[,c('id','Slope')],rep(NA,nrow(ten.day)),rep(0,nrow(ten.day)))
colnames(pred)	<-	c('id.nl2','nl2','id.slope','Slope','new.idnl2','new.nl2')

for (i in 1:nrow(pred)) #Combine DIP rate and 10d data into one matrix
{
	pred[grep(pred[i,1],pred[,3]),c('new.idnl2','new.nl2')] <- pred[i,c('id.nl2','nl2')]
	#pred[pred[,3] %in% pred[i,1],c('new.rank','new.slope')]	<-	
}

dev.new(width=3,height=4)
fn6	<-	paste(write.dir,"/Fig 4 outcome correlation.pdf",sep="")
#pdf(fn6,width=3, height=4)
par(cex=1.1)
plot(pred$Slope,pred$new.nl2,xlab="Prolif rate",ylab="Fold change colony size",
	xlim=c(-.05,.25))
pred.lm	<-	coef(lm(pred$new.nl2~pred$Slope))
abline(pred.lm[1],pred.lm[2])
text(x=.05,y=2.3,paste("R = ",as.numeric(round(cor.test(pred$new.nl2,pred$Slope)$estimate,2))))
#dev.off()

##################################################################
#Code to generate waterfall plots of the data

WF.plot	<-	function(data,yl,titles)
{
	data <- data[order(-data)]
	#print(data)
	plot(0, 0, ylim=yl, xlim=c(0, length(data)),main=titles,
		typ='n', xlab='', ylab='Growth Rate', xaxt='n')
    
	for(i in 1:length(data)) 
		{lines(c(i, i), c(0, data[i]), lwd=150/length(data), col=ifelse(i%%2==0,'gray','black'))}
	abline(h=0, col='red')
	text(x=length(data)*0.8,y=yl[2]*0.9,paste("n = ",length(data),sep=""))
}

dev.new(width=6,height=4)
fn5	<-	paste(write.dir,"/Fig 4 WF plots.pdf",sep="")
#pdf(fn5,width=6, height=4)
par(mfrow=c(1,2))#, mar=c(2,2,2,2)+0.1,font.lab=2)
WF.plot(D.slope$Slope,c(-0.5,1.5),"DMSO")
WF.plot(DCHX.slope$Slope,c(-0.5,1.5),"500 nM CHX")
#dev.off()

##################################################################
#Code to generate a correlation plot of rate and 10-day change in colony size (outcome): Supplementary Figure

ten.day	<-	subset(D.CHX500,Time.day==10)
pred	<-	cbind(ten.day[,c('id','l2')],DCHX.slope[,c('id','Slope')],rep(NA,nrow(ten.day)),rep(0,nrow(ten.day)))
colnames(pred)	<-	c('id.l2','l2','id.slope','Slope','new.idl2','new.l2')

for (i in 1:nrow(pred))
{
	pred[grep(pred[i,1],pred[,3]),c('new.idl2','new.l2')] <- pred[i,c('id.l2','l2')]
	#pred[pred[,3] %in% pred[i,1],c('new.rank','new.slope')]	<-	
}

dev.new(width=3,height=4)
fn6	<-	paste(write.dir,"/Fig 4 outcome correlation.pdf",sep="")
#pdf(fn6,width=3, height=4)
par(cex=1.1)
plot(pred$Slope,pred$new.l2,xlab="Prolif rate",ylab="Fold change colony size",
	ylim=c(-.5,2.4),xlim=c(-.05,.25))
pred.lm	<-	coef(lm(pred$new.l2~pred$Slope))
abline(pred.lm[1],pred.lm[2])
text(x=.05,y=2.3,paste("R = ",as.numeric(round(cor.test(pred$new.l2,pred$Slope)$estimate,2))))
#dev.off()

####################################################################
#Code to examine the error in image processing

e	<-	subset(cfp,Manual.ct>1)
e	<-	subset(e,id!="B02_R_E0_CHX500_TRM0_X17A0_Plate1")
e$lMan.ct	<-	log2(e$Manual.ct)
e$lAuto.ct	<-	log2(e$Cell.count)

elm	<-	lm(e$lMan.ct ~ e$lAuto.ct)
dev.new(width=3,height=4)
fn9	<-	paste(write.dir,"/Fig 2 Correlation ImProcessing.pdf",sep="")
#pdf(fn9,width=3, height=4)
par(font.lab=2)
plot(e$lAuto.ct,e$lMan.ct,ylab="Manual counts",xlab="Automatic counts",ylim=c(3.75,7.75),xlim=c(3.75,7.75))
text(x=5.5,y=7.3,paste("adj. R2 = ",round(summary(elm)$adj.r.squared,2),sep=""))
text(x=5,y=7,paste("n = ",nrow(e),sep=""))
abline(coefficients(elm)[1],coefficients(elm)[2])
#dev.off()

dev.new(width=4,height=4)
fn7	<-	paste(write.dir,"/Fig 2 Residuals ImProcessing.pdf",sep="")
#pdf(fn7,width=4, height=4)
plot(resid(elm),ylim=c(-0.4,0.4),xlab="fitted values",ylab="residuals")
#dev.off()

#########################################################################
#Residuals to linear model fits (unused)


x<-subset(D.CHX500,id==IDCHX500[1])
res<-resid(lm(x$Time.day~x$l2))

for (i in 2:nrow(D.CHX500))
{
	x<-subset(D.CHX500,id==IDCHX500[i])
	res<-rbind(res,resid(lm(x$Time.day~x$l2)))
}
res	<-	data.frame(res)

############################################################################
#QQ plot of image processing error

dev.new(width=3,height=4)
fn10	<-	paste(write.dir,"/Fig 2 QQ-plot ImProcessing.pdf",sep="")
#pdf(fn10,width=3, height=4)
par(font.lab=2)
#qqnorm(DCHX.slope$Slope,xlim=c(-2.5,2.5),xlab="Theoretical quantiles",ylab="Sample quantiles")
qqnorm(DCHX.slope$Slope,xlab="Theoretical quantiles",ylab="Sample quantiles")
qqline(DCHX.slope$Slope)
#dev.off()

