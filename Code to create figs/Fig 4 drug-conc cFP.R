
basefile	<-	'~/Dropbox/Shared Vito/cfp paper/Figures/'
read.dir	<-	paste(basefile,'Data used for figures (pulled by R)',sep='')
setwd(read.dir)
write.dir	<-	paste(basefile,'Figure parts')


source('Load cfp data.R')
source('HG model.R')
library(sn)

#Populate new vectors of DIP rates	
D.cfp	<-	subset(cfp.rates, grepl('D_only',ID))$rates

TRM	<-	c('D.TRM500.cfp','D.TRM50.cfp','D.TRM5.cfp')
D.TRM500.cfp<-	subset(cfp.rates, grepl('_D_',ID) & grepl('TRM500_',ID))$rates
D.TRM50.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('TRM50_',ID))$rates
D.TRM5.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('TRM5_',ID))$rates

FSK	<-	c('D.FSK10.cfp','D.FSK1.cfp','D.FSK0.1.cfp')
D.FSK10.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('FSK10_',ID))$rates
D.FSK1.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('FSK1_',ID))$rates
D.FSK0.1.cfp<-	subset(cfp.rates, grepl('_D_',ID) & grepl('FSK0.1_',ID))$rates

CHX	<-	c('D.CHX500.cfp','D.CHX50.cfp','D.CHX5.cfp')
D.CHX500.cfp<-	subset(cfp.rates, grepl('_D_',ID) & grepl('CHX500_',ID))$rates
D.CHX50.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('CHX50_',ID))$rates
D.CHX5.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('CHX5_',ID))$rates

SB	<-	c('D.SB10.cfp','D.SB1.cfp','D.SB0.1.cfp')
D.SB10.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('SB10_',ID))$rates
D.SB1.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('SB1_',ID))$rates
D.SB0.1.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('SB0.1_',ID))$rates

An	<-	c('D.An1.cfp','D.An0.1.cfp')
D.An1.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('An1_',ID))$rates
D.An0.1.cfp	<-	subset(cfp.rates, grepl('_D_',ID) & grepl('An0.1_',ID))$rates


##########################################################################################
# Add in drug concentration dilution data
# fits with data hidden
##########################################################################################

dev.new(width=3,height=4)
#fn1	<-	paste(write.dir,'Fig 5 p38 pred and valdn.pdf',sep='/')
#pdf(file=fn.p38,width=6,height=4)
#par(font.lab=2,mfrow=c(1,2))
compare.hist(ref='D',combo=c('CHX500','CHX50','CHX5'),my.cols=c('green','red','blue'),
	my.title='CHX',my.xlim=c(-0.05,0.05),my.ylim=c(0,140))

dev.new(width=3,height=4)
#fn1	<-	paste(write.dir,'Fig 5 p38 pred and valdn.pdf',sep='/')
#pdf(file=fn.p38,width=6,height=4)
#par(font.lab=2,mfrow=c(1,2))
compare.hist(ref='D',combo=c('FSK10','FSK1','FSK0.1'),my.cols=c('green','red','blue'),
	my.title='FSK',my.xlim=c(-0.05,0.05),my.ylim=c(0,140))

dev.new(width=3,height=4)
#fn1	<-	paste(write.dir,'Fig 5 p38 pred and valdn.pdf',sep='/')
#pdf(file=fn.p38,width=6,height=4)
#par(font.lab=2,mfrow=c(1,2))
compare.hist(ref='D',combo=c('TRM500','TRM50','TRM5'),my.cols=c('green','red','blue'),
	my.title='TRM',my.xlim=c(-0.05,0.05),my.ylim=c(0,140))

dev.new(width=3,height=4)
#fn1	<-	paste(write.dir,'Fig 5 p38 pred and valdn.pdf',sep='/')
#pdf(file=fn.p38,width=6,height=4)
#par(font.lab=2,mfrow=c(1,2))
compare.hist(ref='D',combo=c('SB10','SB1','SB0.1'),my.cols=c('green','red','blue'),
	my.title='SB',my.xlim=c(-0.05,0.05),my.ylim=c(0,140))

dev.new(width=3,height=4)
#fn1	<-	paste(write.dir,'Fig 5 p38 pred and valdn.pdf',sep='/')
#pdf(file=fn.p38,width=6,height=4)
#par(font.lab=2,mfrow=c(1,2))
compare.hist(ref='D',combo=c('An1','An0.1'),my.cols=c('green','red'),
	my.title='An',my.xlim=c(-0.05,0.05),my.ylim=c(0,140))



##########################################################################################
# Example plot to demonstrate how the curve fits the data
##########################################################################################

dev.new(width=3,height=4)
plot.HG.hist(D.CHX500.cfp,x.limit=c(-0.01,0.05),hist.col='grey',new.plot=T)
plot.HG.hist(D.cfp,x.limit=c(-0.01,0.05),hist.col='white')
legend('topright',c('DMSO','CHX500'),fill=c('white','grey'),cex=0.8,bty='n')


##########################################################################################
# Full plots with data included
##########################################################################################

plot.diln	<-	function(d,my.col,my.ylim=c(0,140),my.linecol)
{
	for (i in 1:length(d))
	{
		print(d[i])
		x	<-	eval(parse(text=d[i]))
		#plot.HG.hist(ref,x.limit=c(-0.05,0.05),hist.col=my.col[1],new.plot=T,show.hist=F,y.limit=my.ylim)
		plot.HG.hist(x,x.limit=c(-0.06,0.06),hist.col=my.col[i],new.plot=T,y.limit=my.ylim,line.col=my.linecol[i])
		legend('topleft',c('ctrl',
			paste(d[i],' n=',length(x))),
			fill=my.col,bty='n',cex=0.9)
		text(x=-0.035,y=100,paste(expression(mu),	format(coef(selm(x~1))['mean'],scientific=T,digits=2),sep='='))
		text(x=-0.035,y=92,paste(expression(sigma),	format(coef(selm(x~1))['s.d.'],scientific=T,digits=2),sep='='))
		text(x=-0.035,y=84,paste(expression(gamma),	format(coef(selm(x~1))[3],scientific=T,digits=2),sep='='))
		}
}


dev.new(width=9,height=4)
par(mfrow=c(1,4))
plot.diln(d=append(CHX,'D.cfp'),my.col=c(rep('grey',3),'white'),my.linecol=c('green','red','blue','black'))


###	This is obsolete now...
dev.new(width=9,height=4)
par(mfrow=c(1,3))
plot.diln(D.cfp,TRM,c('white','grey'))

dev.new(width=9,height=4)
par(mfrow=c(1,3))
plot.diln(D.cfp,FSK,c('white','grey'))

dev.new(width=9,height=4)
par(mfrow=c(1,3))
plot.diln(D.cfp,CHX,c('white','grey'))


dev.new(width=9,height=4)
par(mfrow=c(1,3))
plot.diln(D.cfp,SB,c('white','grey'))

