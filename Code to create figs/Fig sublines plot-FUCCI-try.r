
###For processing of day3-6 subline erlotinib data
setwd('C:/Documents and Settings/BD/My Documents/Dropbox/lab/R-code/BD Population level data');
#setwd('C:/Users/Peter/Dropbox/Lab/R-code/BD Population level data'); for use on laptop 
#PC9.DS <- read.csv('Summary PC9 DS 3-6day erlotinib.csv');
PC9.3day <- read.csv('PC9 sublines full 72hour data_b.csv',row.names=NULL);
head(PC9.3day);
###creates a new column with the log2 data
PC9.3day$log2 <- log2(PC9.3day$Cell.count);
#PC9.3day$ID <- paste(paste(PC9.DS$Subline,PC9.DS$Expt,sep="_"),PC9.DS$Well,sep="_")
###Creates unique ID for each well
PC9.3day$ID <- paste(paste(paste(PC9.3day$Subline,PC9.3day$Expt,sep="_"),
			PC9.3day$Well,sep="_"),PC9.3day$Condition,sep="_")									
###Compiles list of all unique events
well	<-	as.character(unique(PC9.3day$Well))
expt	<-	as.character(unique(PC9.3day$Expt))
subl	<-	as.character(unique(PC9.3day$Subline))
condn	<-	as.character(unique(PC9.3day$Condition))
ID	<-	as.character(unique(PC9.3day$ID))

PC9.3day$Time_h <- floor(PC9.3day$Time_h)

###Normalize data to start at 0 from a log-scale
PC9.3day$norm	<-	rep(0)
norm		<-	numeric()
for (i in ID)
	norm	<-	append(norm, 	subset(PC9.3day, ID == i)$log2 - subset(PC9.3day, ID == i)$log2[1])
PC9.3day$norm	<-	norm
PC9.3day$FUCCI.pos	<- PC9.3day$FUCCI.count/PC9.3day$Cell.count*100
PC9.3day$norm.FUCCI	<-	PC9.3day$norm*PC9.3day$FUCCI.pos/100
###Creates new data files for each subline and condition, e.g., DS1.DMSO, DS4.Erlotinib, etc.
for (i in subl)
	for (j in condn)
		assign(sprintf(paste(i, j, sep = ".")), 
			subset(PC9.3day, PC9.3day$Subline == i & PC9.3day$Condition == j))

### DMSO of all 7 sublines and BR1
#Fills empty data points in the DMSO data
for (j in 1:length(subl))
{
	curr.ds = eval(parse(text=sprintf(paste(subl[j], "DMSO", sep = ".")))) 
	len = length(curr.ds[,1])
	s=0
	for(k in 1:len)
		{
		rownames(curr.ds)[k] <- as.character(k)
		}
	co=0
	for (i in 2:len) #code to find where time increases incorrectly
		if (curr.ds$Time_h[i] > (curr.ds$Time_h[i-1] + 1))
			{#print(i)
			co=co+1}
			print(co)
	for (i in 2:(len+co)) #code to scan data to find where time increases incorrectly
		if (curr.ds$Time_h[i] > (curr.ds$Time_h[i-1] + 1))
			#print(i)
			{
			curr.ds[len+1+s,]<-rep(0,12)
			curr.ds[len+1+s,c(4,5,6,8,9,10,11,12)] <- rep(0) 			#extend data frame by one to allow cut and paste
			curr.ds[len+1+s,c(1,2,3,7)] <- curr.ds[len+s,c(1,2,3,7)] 	#extend data frame by one to allow cut and paste
			row.names(curr.ds)[len+1+s] <- as.character(as.numeric(row.names(curr.ds[len +s,]))+1) #code to advance row names
			curr.ds[(i+1):(len+1+s),] <- curr.ds[(i):(len+s),] 			#code to shift all data ahead 1 hour
			curr.ds[i,c(4,5,6,8,10,11,12)] <-rep(999)
			s=s+1
			}
	#Add code to replace "999" with averaged data from timestep immediately before and after empty data cells
	for (i in 2:(len+co))
		if(curr.ds$Time_h[i] == 999)
			{
			curr.ds[i,c(4,5,6,8,10,11,12)] <- (curr.ds[i-1,c(4,5,6,8,10,11,12)] + curr.ds[i+1,c(4,5,6,8,10,11,12)]) /2
			}
	assign(paste0(subl[j],".DMSO"),curr.ds)
}


DSF	<-	data.frame(Time.h=0:71,DS3.mean=rep(0,72),DS5.mean=rep(0,72),DS3.err=rep(0,72),DS5.err=rep(0,72),DS3.n=rep(0,72),DS5.n=rep(0,72))

for (i in 0:71)
{
	DS3	<-	subset(DS3.DMSO,Time_h==i)$FUCCI.pos
	DS5	<-	subset(DS5.DMSO,Time_h==i)$FUCCI.pos
	DSF$DS3.mean[i+1]	<-	mean(DS3)
	DSF$DS5.mean[i+1]	<-	mean(DS5)
	
	DSF$DS3.95[i+1]	<-	t.test(DS3)$conf.int[2]-mean(DS3)
	DSF$DS5.95[i+1]	<-	t.test(DS5)$conf.int[2]-mean(DS5)
	
	DSF$DS3.SE[i+1]	<-	sd(DS3)/sqrt(length(DS3))
	DSF$DS5.SE[i+1]	<-	sd(DS5)/sqrt(length(DS5))
}

dev.new(width=3,height=4)
plotCI(
x=	DSF$Time.h,
y=	DSF$DS3.mean,
#uiw=DSF$DS3.SE,
ylim=c(0,70)
)

dev.new(width=3,height=4)
plotCI(
x=	DSF$Time.h,
y=	DSF$DS5.mean,
#uiw=DSF$DS5.SE,
ylim=c(0,70)
)

###	Copy and paste code to run plot.72
#multiple sublines example: plot.72(c("DS1.E","BR1.E"),c("norm"),c(0,3),"SE")
#multiple data types example: plot.72("sp.D",c("norm","norm.FUCCI"),c(0,3),"SE")

plot.72	<- function(sublines,data.type,yl,err,color)
{
len	<-	length(sublines)
#dev.new()
par(font.lab=2)
for (i in 1:len)
{
	curr.ds	<-	eval(parse(text=sprintf(paste(sublines[i], "Corr", sep = "."))))
	#con	<-	ifelse(grepl(".E",sublines[i]),"Erlotinib","DMSO")
	pop1	<-	numeric()
	error1	<-	numeric()
	if (length(data.type)>1)
	{pop2	<-	numeric()
	error2	<-	numeric()}
	for (j in 0:71)
	{	r1	<- eval(parse(text=paste("curr.ds$",data.type[1],"[curr.ds$Time_h==","j]",sep="")))
		r1	<-	r1[!is.na(r1)]
		if(length(data.type)>1)
			{r2 <- eval(parse(text=paste("curr.ds$",data.type[2],"[curr.ds$Time_h==","j]",sep="")))
			r2	<-	r2[!is.na(r2)]
			pop2	<-	append(pop2,mean(r2))
			error2	<-	append(error2,ifelse(err=="SE",
						sd(r2)/sqrt(length(r2)),
						(t.test(r2)[[5]][[1]]-t.test(r2)[[4]][[1]])/2))
			}
		pop1	<-	append(pop1,mean(r1))
		error1	<-	append(error1,ifelse(err=="SE",
						sd(r1)/sqrt(length(r1)),
						(t.test(r1)[[5]][[1]]-t.test(r1)[[4]][[1]])/2))
		
	}
	plotCI(
	x	=	c(0:71),
	y	=	pop1,
	uiw	=	error1,
	add	=	ifelse(i==1,FALSE,TRUE),
	ylim	=	yl,
	gap=0,
	ylab="doublings",
	main=ifelse(length(data.type)>1,sublines[i],paste(substring(sublines[i],1,3),data.type))
	#,col="red"
	,col=color[1]
	)
	ifelse(length(data.type)>1,
	plotCI(
	x	=	c(0:71),
	y	=	pop2,
	uiw	=	error2,
	add	=	TRUE,
	ylim	=	yl,
	gap=0,
	ylab="doublings",
	main=ifelse(length(data.type)>1,sublines[i],paste(substring(sublines[i],1,3),data.type))
	,col=color[2]
	),pop2<-numeric()
	)
}
}

ifelse(length(data.type)>1,sublines[i],substring(sublines[i],3)

#Code to plot total cell doublings and FUCCI-only doublings of the same subline and condition per graph
#Plots a subset of the data, with 95% confidence intervals 
plot.subl	<-	c("BR1.D","BR1.E","DS1.D","DS1.E","DS3.D","DS3.E","DS4.D","DS4.E","DS7.D","DS7.E")
dev.new(width=5,height=10)
par(mfrow=c(5,2),mai=c(0.4,0.7,0.1,0.2))
for (i in plot.subl){
	plot.72(i,c("norm","norm.FUCCI"),c(0,3),"S",c("black","dark green"))
	}

#Code to plot the erlotinib-treated and the DMSO treated data of the same data type (either norm cell number or norm-FUCCI) on each graph
#Plots a subset of the data, with 95% confidence intervals
plot.subl	<-	c("BR1","DS1","DS3","DS4","DS7")
dev.new(width=5,height=10)
par(mfrow=c(5,2),mai=c(0.4,0.7,0.1,0.2))
for (i in plot.subl){
	x	<-	c(paste(i,".D",sep=""),paste(i,".E",sep=""))
	plot.72(x,c("norm"),c(0,3),"S",c("black","dark green"))
	plot.72(x,c("norm.FUCCI"),c(0,3),"S",c("black","dark green"))
	}
	
plot.72(c("DS1.D","DS1.E"),c("norm.FUCCI"),c(0,3),"SE")
###############above code only works for "SE"	
	
#multiple sublines example: plot.72(c("DS1.E","BR1.E"),c("norm"),c(0,3),"SE")