basefile	<-	'~/Dropbox/Shared Vito/cFP paper/Figures/'#Code to automatically add PC compatibility
basefile2	<-	gsub('\\\\','/',path.expand(basefile))
xx			<-	try(setwd(basefile2))
if(class(xx)=="try-error"){basefile2 <- paste('C:/Users/Peter',substring(basefile,2),sep="")}
read.dir 	<- paste(basefile2,"Data used for figures (pulled by R)",sep="")
write.dir	<- paste(basefile2,"Figure parts",sep="")
setwd(read.dir)

##################################################################################################
# Read in all data

d	<-	read.csv('PC9_DMSO_FUCCI_Col.csv')
d$id	<-	paste(d$Expt,d$Well,d$Colony,sep="_")
d$FUCCI.pos	<-	round(d$FUCCI/d$Cell.count*100,2)
d$l2	<-	log2(d$Cell.count)
d$nl2	<-	rep(0)
for (i in unique(d$id))
	{x	<-	subset(d,id==i)
	d[grep(i,d$id),'nl2']	<-	x$l2-x$l2[1]
	}

C06L	<-	subset(d,id=="X20120608_C06_L")
D07A	<-	subset(d,id=="X20120410_D07_A")


dev.new(width=3,height=4)
par(font.lab=2)
fn1	<-	paste(write.dir,"/Fig 4 DMSO Colony FUCCI Dynamics.pdf",sep="")
#pdf(fn1,width=3, height=4)
plot(C06L$Time,C06L$FUCCI.pos,ylim=c(0,70),type='o',ylab='Percent S/G2/M',
	xlab='days',pch=16,col='royalblue')
points(D07A$Time,D07A$FUCCI.pos,type='o',pch=21,col='darkorange3')
#dev.off()

dev.new(width=3,height=4)
par(font.lab=2)
fn2	<-	paste(write.dir,"/Fig 4 DMSO Colony Growth Dynamics.pdf",sep="")
#pdf(fn2,width=3, height=4)
plot(C06L$Time,C06L$nl2,ylim=c(0,3),xlim=c(0,3),
	ylab='Colony doublings',xlab='days',pch=16,col='royalblue')
abline(coefficients(lm(nl2 ~ Time,data=C06L)),col='royalblue')
points(D07A$Time,D07A$nl2,pch=21,col='darkorange3')
abline(coefficients(lm(nl2 ~ Time,data=D07A)),col='darkorange3')
#dev.off()


#source('ColonyAssayData.R')
#source('estimate-slopes.R')
#source('estimate-slopes pre-post.R')
#D		<-	subset(slopes,(Condition=="DMSO"&Cell.line=="PC9.par"))
#D$Slope	<-	D$Slope*24/log(2)
#subset(D,Slope>0.85&Slope<0.95)
#subset(D,Slope>0.5&Slope<0.6)	