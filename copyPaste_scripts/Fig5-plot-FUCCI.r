source('find-cFP-Folder.r')
setwd(read.dir)

library(gplots)
library(scales)

PC9.3day <- read.csv('PC9 sublines full 72hour data_b.csv',row.names=NULL);
PC9.3day$log2 <- log2(PC9.3day$Cell.count)		#creates a new column with the log2 data

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

fn	=	paste0(write.dir,'Fig 5 Sublines DMSO per-FUCCI.pdf')
dev.new(width=5,height=4)
#pdf(file=fn,width=5,height=4)
par(mfrow=c(1,2),font.lab=2)

plotCI(
x=	DSF$Time.h,
y=	DSF$DS3.mean,
uiw=DSF$DS3.SE,
#uiw=DSF$DS3.95,
ylim=c(0,70),
ylab= 'Percent S/G2/M',
xlab= 'Time (h)',
main='DS3 DMSO',
col=alpha('black',0.35),
gap=0,
sfrac=0.001,
pch=26
)
lines(x=	DSF$Time.h,y=	DSF$DS3.mean,lwd=2)

plotCI(
x=	DSF$Time.h,
y=	DSF$DS5.mean,
uiw=DSF$DS5.SE,
#uiw=DSF$DS5.95,
ylim=c(0,70),
ylab= 'Percent S/G2/M',
xlab= 'Time (h)',
main='DS5 DMSO',
col=alpha('black',0.35),
gap=0,
sfrac=0.001,
pch=26
)
lines(x=	DSF$Time.h,y=	DSF$DS5.mean,lwd=2)
#dev.off()