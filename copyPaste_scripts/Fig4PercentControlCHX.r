source('find-cFP-Folder.r')
library(gplots)

read.dir<-	paste0(base.dir,'Data used for figures (pulled by R)/')
d.WW	<-	read.csv(paste0(read.dir,'CHX whole well processing.csv'))


s		<-	1.075	#Scaling factor to start at 100%
mean.DMSO	<-	numeric()
for (i in 0:10)	mean.DMSO <- append(mean.DMSO,mean(subset(d.WW,Time.day==i)$Count))	
				
perc	<-	numeric()
err		<-	numeric()
for (i in 0:10)
{
	x		<-	(subset(d.WW,CHX==500 & Time.day==i)$Count/s)/mean.DMSO[i+1]*100
	d.CHX	<-	mean(x)
	perc	<-	append(perc,d.CHX)
	err		<-	append(err,sd(x)/sqrt(length(x)))	#This should be corrected if used; smaller error for smaller values is not right.
}

dev.new(width=3,height=4)
#pdf(file=paste0(write.dir,'Fig 4 per-control CHX.pdf'),width=3,height=4)
par(font.lab=2)
plot(perc, xlab='Time in drug (d)', ylab='%Control',type='o',ylim=c(0,110),main='PC9 CHX500')
#dev.off()