basefile	<-	'~/Dropbox/Lab/cFP paper/Figures/'
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
source('cFP combo estimate-slopes(inprogress).r')

##################################################################
#Code to draw boxplots and show representative colony dynamics

D.CHX500	<- subset(cfp, Erl == 0 & CHX == 500)
D.ctrl		<- subset(cfp, Erl == 0 & CHX == 0 & FSK == 0 & TRM == 0 & X17A == 0)
dev.new(width=3, height=4)
fn.boxplot	<-	paste(write.dir,"/Fig 3 Colony dynamics boxplots.pdf",sep="")
#pdf(fn.boxplot,width=3, height=4)
par(font.lab=2)
boxplot(l2 ~ Time.day, data=D.CHX500[D.CHX500$Time.day>2,], notch=TRUE, xlab="Time in drug(d)", ylab="Colony Doublings", 
		border='grey',ylim=c(-1,2))
n	<- length(unique(D.CHX500$id))
text(y=2,x=6,paste("n = ",n,sep=""))

#lo<-unique(cfp[cfp$l2<0.7&cfp$CHX==500&cfp$Time.day==10&cfp$Erl==0,]$id)
#hi<-unique(cfp[cfp$l2>1.7&cfp$CHX==500&cfp$Erl==0,]$id)
#lo	<-	"B03_B_E0_CHX500_TRM0_X17A0_Plate1"
lo	<-	"B04_J_E0_CHX500_TRM0_X17A0_Plate1"
hi	<-	"B03_A_E0_CHX500_TRM0_X17A0_Plate1"
lines(1:8, subset(cfp,id==lo)$l2[4:11], col='blue', lwd=2)
lines(1:8, subset(cfp,id==hi)$l2[4:11], col='red', lwd=2)
mean.ctrl	<-	rep(0,4)
for (i in 1:4) {mean.ctrl[i]<-mean(subset(D.ctrl,Time.day==i-1)$l2)}
mean.ctrl	<-	mean.ctrl-mean.ctrl[1]
lines(1:4, mean.ctrl,col='black',lty=2,lwd=2)
#dev.off()

##################################################################
#Code to draw linear model fits and boxplots

DCHX.slope		<- subset(slopes,Erl==0&CHX==500)


dev.new(width=3,height=4)
fn4	<-	paste(write.dir,"/Fig 3 linear model colony dynamics.pdf",sep="")
#pdf(fn4,width=3, height=4)
par(font.lab=2)
plot(	4:11, 4:11, type='n', 
		xlab="Time in drug (d)", ylab="Colony Doublings",
		xlim=c(3,10),ylim=c(-1,2),xaxt='n')
for (i in 1:nrow(DCHX.slope))
{
	segments(3,0,10,7*DCHX.slope$Slope[i],col='grey')
}
lm.lo	<-	coef(lm(subset(cfp,id==lo&Time.day>2)$l2~subset(cfp,id==lo&Time.day>2)$Time.day))[2]
lm.hi	<-	coef(lm(subset(cfp,id==hi&Time.day>2)$l2~subset(cfp,id==hi&Time.day>2)$Time.day))[2]
sn		<-	s[4:11]-s[4]
lm.mean	<-	coef(lm.mean<-lm(sn~c(3:10)))[2]
segments(3,0,10,7*lm.hi,col='red',lwd=2)
segments(3,0,10,7*lm.lo,col='blue',lwd=2)
segments(3,0,10,7*lm.mean,col='black',lwd=2)
segments(3,0,6, 3*mean(D.slope$Slope),col='black',lty=2,lwd=2)
#dev.off()

est.col	<-	matrix(NA,nrow=0,ncol=2)
colnames(est.col)	<-	c("Col.size","Time.day")
for (i in c(1:8)){
	x1	<-	DCHX.slope$Slope*(i-1)		
	#x2	<-	rep(i+2,n)
	x2	<-	rep(i+2,n)
	est.col		<-	rbind(est.col,cbind(x1,x2))
}
#boxplot(Col.size ~ Time.day,data=est.col,add=TRUE)

dev.new(width=3,height=4)
par(font.lab=2)
boxplot(Col.size ~ Time.day,data=est.col,
	xlab="Time in drug (d)", ylab="Colony Doublings",
		ylim=c(-1,2))

r	<-	2 #registration error of boxplots		
for (i in 1:nrow(DCHX.slope))
{
	segments(3-r,0,10-r,7*DCHX.slope$Slope[i],col='grey')
}
lm.lo	<-	coef(lm(subset(cfp,id==lo&Time.day>2)$l2~subset(cfp,id==lo&Time.day>2)$Time.day))[2]
lm.hi	<-	coef(lm(subset(cfp,id==hi&Time.day>2)$l2~subset(cfp,id==hi&Time.day>2)$Time.day))[2]
sn		<-	s[4:11]-s[4]
lm.mean	<-	coef(lm.mean<-lm(sn~c(3:10-r)))[2]
segments(3-r,0,10-r,7*lm.hi,col='red',lwd=2)
segments(3-r,0,10-r,7*lm.lo,col='blue',lwd=2)
segments(3-r,0,10-r,7*lm.mean,col='black',lwd=2)
segments(3-r,0,6-r, 3*mean(D.slope$Slope),col='black',lty=2,lwd=2)
#dev.off()



##### re-try plot order


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
	segments(3-r,0,10-r,7*DCHX.slope$Slope[i],col='grey')
}
lm.lo	<-	coef(lm(subset(cfp,id==lo&Time.day>2)$l2~subset(cfp,id==lo&Time.day>2)$Time.day))[2]
lm.hi	<-	coef(lm(subset(cfp,id==hi&Time.day>2)$l2~subset(cfp,id==hi&Time.day>2)$Time.day))[2]
sn		<-	s[4:11]-s[4]
lm.mean	<-	coef(lm.mean<-lm(sn~c(3:10)))[2]
segments(3-r,0,10-r,7*lm.hi,col='red',lwd=2)
segments(3-r,0,10-r,7*lm.lo,col='blue',lwd=2)
segments(3-r,0,10-r,7*lm.mean,col='black',lwd=2)
segments(3-r,0,6-r, 3*mean(D.slope$Slope),col='black',lty=2,lwd=2)
#dev.off()

est.col	<-	matrix(NA,nrow=0,ncol=2)
colnames(est.col)	<-	c("Col.size","Time.day")
for (i in c(1:8)){
	x1	<-	DCHX.slope$Slope*(i-1)		
	#x2	<-	rep(i+2,n)
	x2	<-	rep(i-1,n)
	est.col		<-	rbind(est.col,cbind(x1,x2))
}

boxplot(Col.size ~ Time.day,data=est.col,add=TRUE,xaxt='n')
axis(1,at=1:8,labels=c(3:10))


for (i in 1:nrow(DCHX.slope))
{
	segments(3-r,0,10-r,7*DCHX.slope$Slope[i],col='grey')
}
lm.lo	<-	coef(lm(subset(cfp,id==lo&Time.day>2)$l2~subset(cfp,id==lo&Time.day>2)$Time.day))[2]
lm.hi	<-	coef(lm(subset(cfp,id==hi&Time.day>2)$l2~subset(cfp,id==hi&Time.day>2)$Time.day))[2]
sn		<-	s[4:11]-s[4]
lm.mean	<-	coef(lm.mean<-lm(sn~c(3:10-r)))[2]
segments(3-r,0,10-r,7*lm.hi,col='red',lwd=2)
segments(3-r,0,10-r,7*lm.lo,col='blue',lwd=2)
segments(3-r,0,10-r,7*lm.mean,col='black',lwd=2)
segments(3-r,0,6-r, 3*mean(D.slope$Slope),col='black',lty=2,lwd=2)
#dev.off()















