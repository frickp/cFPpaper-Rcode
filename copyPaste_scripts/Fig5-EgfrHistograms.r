importPackages('sn')
importPackages('scales')

#source('HGmodel.r')
source(textConnection(getURL(paste0(mybaseURL,'HGmodel.r'))))

#read.csv(textConnection(getURL(paste0(mybaseURL,'2011-10-27-IF-EgfrSingleCellIntensity.csv'))))

IF	<-	read.csv(textConnection(getURL(paste0(mybaseURL,'2011-10-27-IF-EgfrSingleCellIntensity.csv'))))
importPackages('scales') #for alpha blending

IF.hist	<-	paste(write.dir,"/Fig 5 EGFR IF histograms.pdf",sep="")
dev.new(width=3,height=4)
#pdf(file=IF.hist, width=3,height=4)
par(font.lab=2)
plot.HG.hist(IF$DS3,y.limit=c(0,2e-3),hist.col=alpha('red',0.5),x.limit=c(500,4000),
	new.plot=T,my.bin=100,line.col='red')
plot.HG.hist(IF$DS5[!is.na(IF$DS5)],hist.col=alpha('blue',0.5),x.limit=c(500,4000),
	my.bin=100,line.col='blue')
legend("topright",c("DS3","DS5"),col=c(alpha('red',0.5),alpha('blue',0.5)),pch=15,bty="n")
#dev.off()

#Compute stats for DS3 using exponentially-modified gaussian
library(fracprolif)
m.DS3	<-	emg.mle(IF$DS3) #computes best fits for mu, sig, and lambda
ks.test(IF$DS3,coef(m.DS3)['mu'],coef(m.DS3)['sigma'],coef(m.DS3)['lambda']) #tests the model fit to the data
#curve(demg(x, coef(m.DS3)['mu'], coef(m.DS3)['sigma'], coef(m.DS3)['lambda']), add=TRUE, col='red', lwd=2)
#curve(demg(x, coef(m.DS3)['mu']-1/coef(m.DS3)['lambda'], sd(IF$DS3), coef(m.DS3)['lambda']), add=TRUE, col='red', lwd=2)
curve(demg(x, coef(m.DS3)['mu']-1/coef(m.DS3)['lambda'], 
	sqrt(coef(m.DS3)['sigma']^2-1/(coef(m.DS3)['lambda']^2)), 
	coef(m.DS3)['lambda']), add=TRUE, col='red', lwd=2)

#Compute stats for DS3
m.DS5	<-	emg.mle(IF$DS5[!is.na(IF$DS5)]) #computes best fits for mu, sig, and lambda
ks.test(IF$DS5,coef(m.DS5)['mu'],coef(m.DS5)['sigma'],coef(m.DS5)['lambda']) #tests the model fit to the data
#curve(demg(x, coef(m.DS5)['mu'], coef(m.DS5)['sigma'], coef(m.DS5)['lambda']), add=TRUE, col='blue', lwd=2)

curve(demg(x, coef(m.DS5)['mu']-1/coef(m.DS5)['lambda'], 
	sqrt(coef(m.DS5)['sigma']^2-1/(coef(m.DS5)['lambda']^2)), 
	coef(m.DS5)['lambda']), add=TRUE, col='blue', lwd=2)

#dev.off()
