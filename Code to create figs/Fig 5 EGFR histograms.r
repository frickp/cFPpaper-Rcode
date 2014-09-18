basefile	<-	'~/Dropbox/Shared Vito/cFP paper/Figures/'
#Code to automatically add PC compatibility
basefile2	<-	gsub('\\\\','/',path.expand(basefile))
xx			<-	try(setwd(basefile2))
if(class(xx)=="try-error"){basefile2 <- paste('C:/Users/Peter',substring(basefile,2),sep="")}
read.dir 	<- paste(basefile2,"Data used for figures (pulled by R)",sep="")
write.dir	<- paste(basefile2,"Figure parts",sep="")
setwd(read.dir)

library(sn)
IF	<-	read.csv('2011-10-27 IF EGFR single-cell intensity.csv')
library(scales) #for alpha blending

IF.hist	<-	paste(write.dir,"/Fig 5 EGFR IF histograms.pdf",sep="")
dev.new(width=3,height=4)
#pdf(file=IF.hist, width=3,height=4)
par(font.lab=2)
hist(IF$DS3,freq=F,ylim=c(0,2e-3),col=alpha('red',0.5),breaks=20,xlim=c(500,3500),
	xlab='Fluorescence intensity (A.U.)',main='EGFR intensity')
hist(IF$DS5,add=T,freq=F,col=alpha('blue',0.5),breaks=20)
legend("topright",c("Clone1","Clone2"),col=c(alpha('red',0.5),alpha('blue',0.5)),pch=15,bty="n")

m.DS3	<-	selm(IF$DS3~1)
curve(dsn(x,xi=coef(m.DS3,"DP")['xi'],omega=coef(m.DS3,"DP")['omega'],alpha=coef(m.DS3,"DP")['alpha']),add=T,col='red',lwd=2)
ks.test(IF$DS3[!is.na(IF$DS3)],'psn',xi=coef(m.DS3,"DP")['xi'],omega=coef(m.DS3,"DP")['omega'],alpha=coef(m.DS3,"DP")['alpha'])

m.DS5	<-	selm(IF$DS5~1)
curve(dsn(x,xi=coef(m.DS5,"DP")['xi'],omega=coef(m.DS5,"DP")['omega'],alpha=coef(m.DS5,"DP")['alpha']),add=T,col='blue',lwd=2)
ks.test(IF$DS5[!is.na(IF$DS5)],'psn',xi=coef(m.DS5,"DP")['xi'],omega=coef(m.DS5,"DP")['omega'],alpha=coef(m.DS5,"DP")['alpha'])
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
