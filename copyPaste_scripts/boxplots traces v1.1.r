##################################################################
# Pull all data
##################################################################

importPackages('scales') #for alpha blending
importPackages('gplots') #for error bars: plotCI
importPackages('RCurl') #for error bars: plotCI
source(textConnection(getURL(paste0(mybaseURL,'cFP-norm72h.r'))))
source(textConnection(getURL(paste0(mybaseURL,'cFP-comboEstimateSlopes.r'))))


##################################################################
# Subset the data to highlight one drug

D.CHX500	<- subset(cfp, Erl == 0 & CHX == 500)
IDCHX500	<- unique(D.CHX500$id)

DCHX.slope		<- subset(slopes,Erl==0 & CHX==500)
DCHX.slope50	<- subset(slopes,Erl==0 & CHX==50)
DCHX.slope5		<- subset(slopes,Erl==0 & CHX==5)
D.slope			<- subset(slopes, Erl == 0 & CHX == 0 & FSK == 0 & TRM == 0 & X17A == 0)

##################################################################
#Code to plot the linear model fit r-squared values for CHX-treated colonies: Supplementary Figure

dev.new(width=3,height=4)
par(font.lab=2)
R2	<-	append(D.slope$R2,DCHX.slope$R2)
#R2	<-	append(	append(	append(D.slope$R2, DCHX.slope$R2),DCHX.slope50$R2),DCHX.slope5$R2)
hist(R2,xlab="R2",main="",breaks=20)
text(x=0,y=250,paste("n = ",length(R2),sep=""))


##################################################################
#Code to generate a correlation plot of rate and 10-day change in colony size (outcome): Supplementary Figure

ten.day	<-	subset(D.CHX500,Time.day==10)
pred	<-	cbind(ten.day[,c('id','nl2')],DCHX.slope[,c('id','Slope')],rep(NA,nrow(ten.day)),rep(0,nrow(ten.day)))
colnames(pred)	<-	c('id.nl2','nl2','id.slope','Slope','new.idnl2','new.nl2')

for (i in 1:nrow(pred)) #Combine DIP rate and 10d data into one matrix
{
	pred[grep(pred[i,1],pred[,3]),c('new.idnl2','new.nl2')] <- pred[i,c('id.nl2','nl2')]
}

dev.new(width=3,height=4)
par(cex=1.1)
plot(pred$Slope,pred$new.nl2,xlab="Prolif rate",ylab="Fold change colony size",
	xlim=c(-.05,.25))
pred.lm	<-	coef(lm(pred$new.nl2~pred$Slope))
abline(pred.lm[1],pred.lm[2])
text(x=.05,y=2.1,paste("R = ",as.numeric(round(cor.test(pred$new.nl2,pred$Slope)$estimate,2))))

####################################################################
#Code to examine the error in image processing

e	<-	subset(cfp,Manual.ct>1)
e	<-	subset(e,id!="B02_R_E0_CHX500_TRM0_X17A0_Plate1")
e$lMan.ct	<-	log2(e$Manual.ct)
e$lAuto.ct	<-	log2(e$Cell.count)

elm	<-	lm(e$lMan.ct ~ e$lAuto.ct)
dev.new(width=3,height=4)
par(font.lab=2)
plot(e$lAuto.ct,e$lMan.ct,ylab="Manual counts",xlab="Automatic counts",ylim=c(3.75,7.75),xlim=c(3.75,7.75))
text(x=5.5,y=7.3,paste("adj. R2 = ",round(summary(elm)$adj.r.squared,2),sep=""))
text(x=5,y=7,paste("n = ",nrow(e),sep=""))
abline(coefficients(elm)[1],coefficients(elm)[2])

dev.new(width=4,height=4)
plot(resid(elm),ylim=c(-0.4,0.4),xlab="fitted values",ylab="residuals")


############################################################################
#QQ plot of image processing error

dev.new(width=3,height=4)
fn10	<-	paste(write.dir,"/Fig 2 QQ-plot ImProcessing.pdf",sep="")
par(font.lab=2)
qqnorm(DCHX.slope$Slope,xlab="Theoretical quantiles",ylab="Sample quantiles")
qqline(DCHX.slope$Slope)


