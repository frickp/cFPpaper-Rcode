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

#####################################################################
#Plot average dynamics in 500 nM CHX versus control

ctrl.mean	<-	aggregate(D.ctrl$nl2,by=list(D.ctrl$Time.day),FUN=mean)
colnames(ctrl.mean)	<-	c("Time.day","nl2")
CHX.mean	<-	aggregate(D.CHX500$nl2,by=list(D.CHX500$Time.day),FUN=mean)
colnames(CHX.mean)	<-	c("Time.day","nl2")


dev.new(width=3,height=4)
fn2	<-	paste(write.dir,"/Fig 3 population average dynamics.pdf",sep="")
#pdf(fn2,width=3, height=4)
par(font.lab=2)
plot(CHX.mean$Time.day,CHX.mean$nl2,type='o',ylim=c(0,2),xlab="Time in drug (d)",ylab="Population doublings",lwd=2,pch=26)
lines(ctrl.mean$Time.day,ctrl.mean$nl2,lty=2,lwd=2)
rect(3,-0.2,10.3,2.2,col=alpha('lightblue',0.5),border=NA)
text(x=5,y=0.1,"Linear")
legend("topright",c("Ctrl","Chx"),lty=c(1,2),lwd=c(2,2),bty="n")
#dev.off()

##################################################################
#Code to plot the normalized standard error in colony size over time

err	<-	aggregate(D.CHX500$Cell.count, by=list(D.CHX500$Time.day),
            FUN=function(x)
               {
					sd(x)/sqrt(length(x))
				})
dev.new(width=3,height=4)
fn3	<-	paste(write.dir,"/Fig 3 Increasing SE over time.pdf",sep="")
#pdf(fn3,width=3, height=4)
par(font.lab=2)
#plot(c(0:10),err$x,ylim=c(0,6),ylab="Colony size (Std Error)",xlab="days in 500nM CHX")
plot(c(0:10),err$x/err$x[1],ylim=c(1,3),ylab="Colony size std. error",xlab="Time in drug (d)")
#dev.off()

##################################################################
#Code to plot the linear model fit r-squared values for CHX-treated colonies

DCHX.slope		<- subset(slopes,Erl==0 & CHX==500)
DCHX.slope50	<- subset(slopes,Erl==0 & CHX==50)
DCHX.slope5		<- subset(slopes,Erl==0 & CHX==5)
D.slope			<- subset(slopes, Erl == 0 & CHX == 0 & FSK == 0 & TRM == 0 & X17A == 0)

dev.new(width=3,height=4)
fn8	<-	paste(write.dir,"/Fig 4 R2 linear model coefficients.pdf",sep="")
#pdf(fn8,width=3, height=4)
par(font.lab=2)
R2	<-	append(D.slope$R2,DCHX.slope$R2)
#R2	<-	append(	append(	append(D.slope$R2, DCHX.slope$R2),DCHX.slope50$R2),DCHX.slope5$R2)
hist(R2,xlab="R2",main="",breaks=20)
text(x=0,y=250,paste("n = ",length(R2),sep=""))
#dev.off()

lR2	<- subset(DCHX.slope,R2<0.6)

x	<-	subset(D.CHX500,id=="B05_K_E0_CHX500_TRM0_X17A0_Plate1")
dev.new(width=3,height=4)
plot(x$Time.day,x$l2,ylim=c(-1,2))
l	<-	lm(x$l2[4:11] ~ x$Time.day[4:11])
lines(x=c(3,11),y=c(x$l2[4],x$l2[4]+coef(l)[[2]]*7))

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
#Code to generate a correlation plot of rate and 10-day change in colony size (outcome)

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

