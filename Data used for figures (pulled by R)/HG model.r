########################################################################################
#HG model definitions
########################################################################################

hgm.log <- function(t, mu, sigma, y0=1)
{
    log(y0) + mu*t + sigma*sigma*t*t/2
}

HGM.skew <- function(t, location, scale, shape)
{
	delta <- shape / sqrt(1+shape*shape)
	log(2) + location*t + t*t*scale*scale/2 + pnorm(scale*delta*t, log.p=TRUE)
}

########################################################################################
# HG model plotting function
########################################################################################
#parameter definitions
#data:		is the raw dataset (numbers only) to be fit
#skewness:	either T or F; for skew-normal or normal fits respectively 
#x.limits: 	defines the range of x-values, ex.: c(0,30)
#y.limits: 	defines the range of y-values, ex.: c(-5,5)
#new.plot:	either T or F; if curve is being added to existing plot, then use F
#cols:		color of the curve

plot.HGM	<-	function(data,skewness,x.limits,y.limits,new.plot,cols){
d.mu		<-	mean(data)
d.sigma		<-	ifelse(is.na(sd(data)),0,sd(data))
if(skewness==T){
d.xi		<-	coef(selm(data~1),"DP")['xi']
d.omega		<-	coef(selm(data~1),"DP")['omega']
d.alpha		<-	coef(selm(data~1),"DP")['alpha']
}
if(new.plot==T){
	plot(	0:x.limits[2], 0:x.limits[2], type='n', 
			xlab="Time in erl (d)", ylab="Population doublings",
			xlim=x.limits,ylim=y.limits)
}
if(skewness==T){
	curve(HGM.skew(x*24,d.xi,d.omega,d.alpha),from=0, to=x.limits[2], lwd=3, 
		add=TRUE,col=cols)
	if(new.plot==T)	{title(expression(bold("HG model (skew-normal)")))}	
}
else{
	curve(	hgm.log(x*24, d.mu, d.sigma)/log(2), 
			from=0, to=x.limits[2], lwd=3, 
			add=TRUE,col=cols)
	if(new.plot==T)	{title(expression(bold("HG model (normal)")))}
}
abline(h=0, col='darkgray', lty=2, lwd=2)
}

########################################################################################
# plotting function for histograms used in HG model
########################################################################################
plot.HG.hist	<-	function(d,x.limit=c(-0.05,0.05),y.limit=c(0,150),skewness=T,line.col='black',hist.col,line.type=1,line.width=2,
						new.plot=F,show.hist=T,my.title=NA,my.border='black',my.bin=0.0025)
{
	if(new.plot==T)	
	{
		plot(seq(x.limit[1],x.limit[2],my.bin),seq(x.limit[1],x.limit[2],my.bin),
		xlim=x.limit,ylim=y.limit,type='n',xlab='DIP rate',ylab='Density',main=my.title)
		}
	if(show.hist==T)	
	{
		hist(d,freq=F,breaks=seq(x.limit[1],x.limit[2],my.bin),ylim=y.limit,
		col=hist.col,border=ifelse(show.hist==T,my.border,alpha(my.border,0)),add=T)
	}
	d.xi	<-	coef(selm(d~1),"DP")['xi']
	d.omega	<-	coef(selm(d~1),"DP")['omega']
	d.alpha	<-	coef(selm(d~1),"DP")['alpha']	
	if(skewness==TRUE){
		curve(dsn(x,d.xi,d.omega,d.alpha), xlim=x.limit,ylim=y.limit,col=line.col, 
			lwd=line.width, lty=line.type,add=T,xlab='DIP rate',ylab='',
			from=range(d)[1]-my.bin*2, to=range(d)[2]+my.bin*2)
			cat(my.title, "p= ", ks.test(d, 'psn', d.xi,d.omega,d.alpha)$p.value,'\n')			
	}	
	else{
		curve(dnorm(x, mean(d), sd(d)), xlim=x.limit,ylim=y.limit,col=line.col, 
			lwd=line.width, lty=line.type,add=T,main=name,xlab='DIP rate',ylab='')
	}
}


########################################################################################
# predictions of TTR and DOR
########################################################################################

HGM.TTR	<-	function(data,skewness){
d.xi	<-	coef(selm(data~1),"DP")['xi']
d.omega	<-	coef(selm(data~1),"DP")['omega']
d.alpha	<-	coef(selm(data~1),"DP")['alpha']
if(skewness==T){
	uniroot(function(t) HGM.skew(t, d.xi,d.omega,d.alpha), lower=0.01, upper=1e4)$root
}
else -2*mean(data)/(sd(data)^2)
}

HGM.DOR	<-	function(data,skewness){
d.xi	<-	coef(selm(data~1),"DP")['xi']
d.omega	<-	coef(selm(data~1),"DP")['omega']
d.alpha	<-	coef(selm(data~1),"DP")['alpha']
if(skewness==T){
	t.min	<-	optimize(f=function(t) HGM.skew(t, d.xi,d.omega,d.alpha), 
		interval=c(1,10000),maximum=F)$minimum
	y.min	<-	HGM.skew(t.min,d.xi,d.omega,d.alpha)
	return(c(t.min,y.min))
}
else{
	t.min	<-	optimize(f=function(t) hgm.log(t, mean(data),sd(data),1), 
		interval=c(1,10000),maximum=F)$minimum
	y.min	<-	hgm.log(t.min,mean(data),sd(data),1)
	return(c(t.min,y.min))
}
}

########################################################################################
#Plotting function for comparing skew-normal hist fits (without the data)
########################################################################################
#ref is the single-agent data
#combo is the second agent added to 'ref'

compare.hist	<- function(ref,combo,my.cols,my.title,my.xlim=c(-0.04,0.02),my.ylim=c(0,100))
{
	plot.HG.hist(d=eval(parse(text=paste(ref,'cfp',sep='.'))), new.plot=T,
		show.hist=F,x.limit=my.xlim,y.limit=my.ylim,line.col='black')
	for (i in 1:length(combo))
	{
		cfp		<- eval(parse(text=paste(ref,combo[i],'cfp',sep='.')))
		plot.HG.hist(d=cfp,show.hist=F,line.col=my.cols[i])
	}
	abline(v=0,col='burlywood',lty=2,lwd=2)
	legend('topleft',paste(ref,append('DMSO',combo),sep='+'),fill=append('black',my.cols),bty='n',cex=0.7)
	title(my.title)
}

########################################################################################
#Plotting function for overlaying the data and the HG model predictions together
########################################################################################

plot.HG.valid	<-	function(ref,combo,my.cols)
{
	#Plot erlotinib alone
	plot.reb(data=eval(parse(text=paste(ref,'reb',sep='.'))), x.limits=c(0,33), 
		y.limits=c(-3,0.5), my.pch=1, my.cex=0.7, new.plot=T, my.col=alpha('black',0.6))
	plot.HGM(d=eval(parse(text=paste(ref,'cfp',sep='.'))), skewness=T,x.limits=c(0,35),
		y.limits=c(-2.5,1),new.plot=F,cols='grey')
	#Plot the combinations
	for (i in 1:length(combo))
	{
		cfp		<- eval(parse(text=paste(ref,combo[i],'cfp',sep='.'))) 
		validn	<- eval(parse(text=paste(ref,combo[i],'reb',sep='.')))
		plot.reb(data=validn, x.limits=c(0,33), y.limits=c(-3,0.5), my.pch=1, my.cex=0.7, 
			my.col=alpha(my.cols[i],0.6))
		plot.HGM(d=cfp,skewness=T,x.limits=c(0,35),y.limits=c(-2.5,1),new.plot=F,cols=my.cols[i])
	}
	legend('topleft',c('data','predictions'),pch=c('o','-'),cex=0.6,bty='n')
}
