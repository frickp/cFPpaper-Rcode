importPackages = function(pkgName)
{
	if(eval(parse(text=paste0('library(',pkgName,',logical.return=T,quietly=T)'))))
	{
		print(paste0(pkgName,': ','found'))
	} else {
		print(paste(pkgName,'not found'))
		if(readline("Type T to install package or anything else to abort")=='T'){
			print(paste('not found',pkgName))
			install.packages(pkgName)
			}
		else {
			print('Package install aborted. Type T to install')
		}
	}
}

importPackages('RCurl')		#for pulling raw data from Github
importPackages('sn')		#skew-normal for histogram fits
importPackages('scales') 	#for alpha blending
importPackages('gplots') 	#for error bars: plotCI


#baseURL = c("https://raw.github.com/--username--/--repo-name--/master/")
mybaseURL='https://raw.githubusercontent.com/frickp/cFPpaper-Rcode/master/copyPaste_rawDataAndPreprocessingScripts/'
source(textConnection(getURL(paste0(mybaseURL,'ColonyAssayData.R'))))

read.csv(textConnection(getURL(paste0(mybaseURL,'BP_SKMEL5.csv'))))

https://rawgit.com/frickp/cFPpaper-Rcode/blob/master/Data%20used%20for%20figures%20(pulled%20by%20R)/BP_SKMEL5.csv




source('find-cFP-Folder.r')
setwd(read.dir)

##################################################################
# Pull all data
##################################################################

source('cFP norm (72h).R')
source('cFP combo estimate-slopes.r')
source('HG model.r')
library(sn)		#skew-normal for histogram fits
library(scales) #for alpha blending
library(gplots) #for error bars: plotCI

##########################################################################################
#Subset the data
##########################################################################################

#Colony dynamics
D.CHX500	<- subset(cfp, Erl == 0 & CHX == 500)
D.ctrl		<- subset(cfp, Erl == 0 & CHX == 0 & FSK == 0 & TRM == 0 & X17A == 0)
IDCHX500	<- unique(D.CHX500$id)

#Colony rates
DCHX.slope		<- subset(slopes,Erl==0 & CHX==500)
D.slope			<- subset(slopes, Erl == 0 & CHX == 0 & FSK == 0 & TRM == 0 & X17A == 0)

##########################################################################################
#Plot average dynamics in 500 ng/ml CHX versus control: Fig. 3a
##########################################################################################
ctrl.mean	<-	aggregate(D.ctrl$nl2,by=list(D.ctrl$Time.day),FUN=mean)
colnames(ctrl.mean)	<-	c("Time.day","nl2")
CHX.mean	<-	cbind(
					aggregate(D.CHX500$nl2,by=list(D.CHX500$Time.day),FUN=mean),
					aggregate(D.CHX500$nl2, by=list(D.CHX500$Time.day),FUN=sd)$x
					)
colnames(CHX.mean)	<-	c("Time.day","nl2","sd")

dev.new(width=3,height=4)
fn1	<-	paste(write.dir,"/Fig 3 mean dynamics + SD.pdf",sep="")
#pdf(fn1,width=3, height=4)
par(font.lab=2)
plot(CHX.mean$Time.day,CHX.mean$nl2,type='o',ylim=c(-0.2,2.2),xlab="Time in drug (d)",ylab="Population doublings",lwd=2,pch=26)
lines(ctrl.mean$Time.day,ctrl.mean$nl2,lty=2,lwd=2)

plotCI(
	x	=	CHX.mean$Time.day,
	y	=	CHX.mean$nl2,
	uiw	=	CHX.mean$sd, 
	pch	=	20,
	gap	=	0,
	add	=	TRUE,
	col	=	alpha('black',0.25)
)
legend("topright",c("Ctrl","Chx"),lty=c(2,1),lwd=c(2,2),bty="n")
#dev.off()
#rect(3,-0.2,10.3,2.2,col=alpha('lightblue',0.5),border=NA)
#text(x=5,y=0.1,"Linear")

##########################################################################################
#Plot traces of individual colonies over colony size boxplots: Fig. 3b
##########################################################################################

lo	<-	"B04_J_E0_CHX500_TRM0_X17A0_Plate1"
hi	<-	"B03_A_E0_CHX500_TRM0_X17A0_Plate1"
down<-	"B05_H_E0_CHX500_TRM0_X17A0_Plate1"
dev.new(width=3,height=4)	
fn2	<-	paste(write.dir,"/Fig 3 boxplot traces.pdf",sep="")
#pdf(fn2,width=3, height=4)
par(cex.lab=1.2)	
boxplot(nl2 ~ Time.day, data=D.CHX500, notch=TRUE, xlab="Time in drug (d)", ylab="Colony doublings", 
		border='grey50',ylim=c(-0.2,2.2))
lines(1:11, subset(cfp,id==lo)$nl2[1:11], col='blue', lwd=2)
lines(1:11, subset(cfp,id==hi)$nl2[1:11], col='red', lwd=2) 
lines(1:11, subset(cfp,id==down)$nl2[1:11], col='green', lwd=2)
#dev.off()
##########################################################################################
#72h normalization of plot traces of individual colonies over colony size boxplots: Fig. 3d
##########################################################################################

dev.new(width=3, height=4)
fn.boxplot	<-	paste(write.dir,"/Fig 3 72h norm colony dynamics boxplots.pdf",sep="")
#pdf(fn.boxplot,width=3, height=4)
par(font.lab=2)
boxplot(nl2.72 ~ Time.day, data=D.CHX500[D.CHX500$Time.day>2,], notch=TRUE, xlab="Time in drug(d)", ylab="Colony Doublings", 
		border='grey20',ylim=c(-1,2))
n	<- length(unique(D.CHX500$id))
text(y=2,x=6,paste("n = ",n,sep=""))

lines(1:8, subset(cfp,id==lo)$nl2.72[4:11], col='blue', lwd=2)
lines(1:8, subset(cfp,id==hi)$nl2.72[4:11], col='red', lwd=2)
lines(1:8, subset(cfp,id==down)$nl2.72[4:11], col='green', lwd=2)
mean.ctrl	<-	aggregate(D.ctrl$nl2,by=list(D.ctrl$Time.day),mean)
colnames(mean.ctrl)	<- c('Time.day','nl2')
lines(1:4, mean.ctrl$nl2,col='black',lty=2,lwd=2)
#dev.off()



##########################################################################################
# Code to display estimated data from linear model fits: Fig. 3e
##########################################################################################
r	<-	2 #registration error of boxplots		
dev.new(width=3,height=4)
fn3	<-	paste(write.dir,"/Fig 3 72h norm linear model colony dynamics.pdf",sep="")
#pdf(fn3,width=3, height=4)
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
lm.down	<-	coef(lm(subset(cfp,id==down&Time.day>2)$nl2.72~subset(cfp,id==down&Time.day>2)$Time.day))[2]
lm.DMSO	<-	coef(lm(nl2 ~ Time.day,data = D.ctrl))[2]

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
segments(3-r,0,10-r,7*lm.down,col='green',lwd=2)
segments(3-r,0,6-r, 3*lm.DMSO,col='black',lty=2,lwd=2)
#dev.off()

##########################################################################################
#Code to plot the linear model fit r-squared values for CHX-treated colonies: Supplementary Figure
##########################################################################################

dev.new(width=3,height=4)
fn4	<-	paste(write.dir,"/Fig 4 R2 linear model coefficients.pdf",sep="")
#pdf(fn4,width=3, height=4)
par(font.lab=2)
R2	<-	append(D.slope$R2,DCHX.slope$R2)
#R2	<-	append(	append(	append(D.slope$R2, DCHX.slope$R2),DCHX.slope50$R2),DCHX.slope5$R2)
hist(R2,xlab="R2",main="",breaks=20)
text(x=0,y=250,paste("n = ",length(R2),sep=""))
#dev.off()

##########################################################################################
#Code to generate a correlation plot of rate and 10-day change in colony size (outcome): Supplementary Figure
##########################################################################################

ten.day	<-	subset(D.CHX500,Time.day==10)
pred	<-	cbind(ten.day[,c('id','nl2')],DCHX.slope[,c('id','Slope')],rep(NA,nrow(ten.day)),rep(0,nrow(ten.day)))
colnames(pred)	<-	c('id.nl2','nl2','id.slope','Slope','new.idnl2','new.nl2')

for (i in 1:nrow(pred)) #Combine DIP rate and 10d data into one matrix
{
	pred[grep(pred[i,1],pred[,3]),c('new.idnl2','new.nl2')] <- pred[i,c('id.nl2','nl2')]
}

dev.new(width=3,height=4)
fn5	<-	paste(write.dir,"/Fig 4 outcome correlation.pdf",sep="")
#pdf(fn5,width=3, height=4)
par(cex=1.1, font.lab=2)
plot(pred$Slope,pred$new.nl2,xlab="Prolif rate",ylab="Fold change colony size",
	xlim=c(-.05,.25),ylim=c(-0.2,2.4))
pred.lm	<-	coef(lm(pred$new.nl2~pred$Slope))
abline(pred.lm[1],pred.lm[2])
text(x=.05,y=2.3,paste("R = ",as.numeric(round(cor.test(pred$new.nl2,pred$Slope)$estimate,2))))
#dev.off()

dev.new(width=3,height=4)
fn.10dh	<-	paste(write.dir,"/Fig 3 10d hist 500CHX.pdf",sep="")
#pdf(fn.10dh,width=3, height=4)
par(font.lab=2)
plot.HG.hist(ten.day$nl2,x.limit=c(-1,3),hist.col=alpha('black',0.3),new.plot=T,my.bin=0.26,y.limit=c(0,1))
#dev.off()
##########################################################################################
# Plot rate histograms
##########################################################################################


dev.new(width=3,height=4)
fn.CHXrates	<-	paste(write.dir,"/Fig 3 CHX rates hist.pdf",sep="")
#pdf(fn.CHXrates,width=3, height=4)
plot.HG.hist(d=DCHX.slope$Slope,new.plot=T,x.limit=c(-0.15,0.3),
	hist.col=alpha('black',0.3),y.limit=c(0,10),my.bin=0.025,skewness=T)
arrows(lm.lo,10,lm.lo,5,length=0.15,lwd=2,col='blue')
arrows(lm.hi,10,lm.hi,5,length=0.15,lwd=2,col='red')
arrows(lm.down,10,lm.down,5,length=0.15,lwd=2,col='green')
#dev.off()

dev.new(width=4,height=4)
fn.hist2	<-	paste(write.dir,"/Fig 3 hist (DMSO).pdf",sep="")
#pdf(fn.hist2,width=4, height=4)
plot.HG.hist(d=D.slope$Slope,new.plot=T,x.limit=c(0.35,1.25),
	hist.col=alpha('black',0.3),y.limit=c(0,6),my.bin=0.03)
#dev.off()