
library(gplots)
source('find-cFP-Folder.r')
setwd(read.dir)

PC9.3day <- read.csv('PC9 sublines full 72hour data.csv',row.names=NULL);
###creates a new column with the log2 data
PC9.3day$log2 <- log2(PC9.3day$Cell.count);
###Creates unique ID for each well
PC9.3day$ID <- paste(paste(paste(PC9.3day$Subline,PC9.3day$Expt,sep="_"),
			PC9.3day$Well,sep="_"),PC9.3day$Condition,sep="_")									
###Compiles list of all unique events
well	<-	as.character(unique(PC9.3day$Well))
expt	<-	as.character(unique(PC9.3day$Expt))
subl	<-	as.character(unique(PC9.3day$Subline))
condn	<-	as.character(unique(PC9.3day$Condition))
ID		<-	as.character(unique(PC9.3day$ID))

PC9.3day$Time_h <- floor(PC9.3day$Time_h)

###Normalize data to start at 0 from a log-scale
PC9.3day$norm	<-	rep(0)
norm		<-	numeric()
for (i in ID)
	norm	<-	append(norm, 	subset(PC9.3day, ID == i)$log2 - subset(PC9.3day, ID == i)$log2[1])
PC9.3day$norm	<-	norm


###Creates new data files for each subline and condition, e.g., DS1.DMSO, DS4.Erlotinib, etc.
for (i in subl)
	for (j in condn)
		assign(sprintf(paste(i, j, sep = ".pre.")), 
			subset(PC9.3day, PC9.3day$Subline == i & PC9.3day$Condition == j))

# Few timepoints were intermittently missing due to slight differences in the image acquisisiton rate 
# This code finds empty entries and replaces them with the average of the two adjacent timepoints
# from the same sample. Data interpolation manually verified.
### DMSO of all 7 sublines
for (j in 1:7)
{
	curr.ds = eval(parse(text=sprintf(paste(subl[j+1], "DMSO", sep = ".pre.")))) 
	len = length(curr.ds[,1])
	s=0
	for(k in 1:len)
		{
		rownames(curr.ds)[k] <- as.character(k)
		}
	co=0
	for (i in 2:len) #code to find where time increases incorrectly
		if (curr.ds$Time_h[i] > (curr.ds$Time_h[i-1] + 1))
			{
			co=co+1
			}
			print(co)
	for (i in 2:(len+co)) #code to scan data to find where time increases incorrectly
		if (curr.ds$Time_h[i] > (curr.ds$Time_h[i-1] + 1))
			{
			curr.ds[len+1+s,]<-rep(0,10)
			curr.ds[len+1+s,c(4,5,6,8,9,10)] <- rep(0) 					#extend data frame by one to allow cut and paste
			curr.ds[len+1+s,c(1,2,3,7)] <- curr.ds[len+s,c(1,2,3,7)] 	#extend data frame by one to allow cut and paste
			row.names(curr.ds)[len+1+s] <- as.character(as.numeric(row.names(curr.ds[len +s,]))+1) 		#code to advance row names
			curr.ds[(i+1):(len+1+s),] <- curr.ds[(i):(len+s),] 			#code to shift all data ahead 1 hour
			curr.ds[i,c(4,5,6,8,10)] <-rep(999)
			s=s+1
			}
	#Add code to replace "999" with averaged data from timestep immediately before and after empty data cells
	for (i in 2:(len+co))
		if(curr.ds$Time_h[i] == 999)
			{
			curr.ds[i,c(4,5,6,8,10)] <- (curr.ds[i-1,c(4,5,6,8,10)] + curr.ds[i+1,c(4,5,6,8,10)]) /2
			}
	assign(paste(paste("DS",j,sep=".pre."),".Corr",sep=".D"),curr.ds)
}


###Creates data frame to input the rates from all conditions (DMSO, pre72, post72) calculated by linear model fits
E.num      <-  length(as.character(unique(DS1.pre.Erlotinib$ID)))   #Number of Erlotinib wells per subline 
D.num      <-  length(as.character(unique(DS1.pre.DMSO$ID)))        #Number of DMSO wells per subline

n <- 13 #Maximum sample names plus data summaries
DS.rates <- data.frame(DS1.D=rep(0,n), DS2.D=rep(0,n), DS3.D=rep(0,n), DS4.D=rep(0,n), DS5.D=rep(0,n), DS6.D=rep(0,n),
				DS7.D=rep(0,n)
			)

#Code that uses 'lm' to calculate the DMSO rates for all wells in each subline
subl	<-	as.character(unique(PC9.3day$Subline))
subl <- subl[2:8] #Keep only DS1-7
for (i in 1:length(subl))
{
	for (j in 1:D.num)
	{
	curr.ds       =  eval(parse(text=sprintf(paste(subl[i], "DMSO", sep = ".pre."))))
	DS.ID         <- unique(curr.ds$ID)
	DS.rates[j,i] <- lm(curr.ds$norm[curr.ds$ID == DS.ID[j] & curr.ds$Time_h < 72] 
				~ curr.ds$Time_h[curr.ds$ID == DS.ID[j] & curr.ds$Time_h < 72])[[1]][[2]]
	}
}

#Code to add mean and 95% confidence error magnitude for each rate
for (i in 1:length(DS.rates))
{
	Col <- DS.rates[,i]
	repl <- Col[Col!=0]									#replicate wells
	DS.rates[11,i] <- mean(repl)						#mean
	DS.rates[12,i] <- t.test(repl)[[4]][[2]]-mean(repl)	#95% confidence
	DS.rates[13,i] <- sd(repl)/sqrt(length(repl))		#standard error
}



#fn1	=	paste(write.dir,"/Fig 4 DMSO subline rates.pdf",sep="")
dev.new(width=2.5,height=4)
#pdf(file=fn1, width=2.5, height=4)
par(font.lab=2)
plotCI(
	barplot(
		as.numeric(DS.rates[11,c(5,3)])*24,
		ylim=c(0,1),
		ylab="Prolif. rate",
		names=c("c5","c3"),
		xlab="days"
	),
	as.numeric(DS.rates[11,c(5,3)])*24,
	uiw	=	as.numeric(DS.rates[13,c(5,3)])*24,
	add	=	TRUE,
	pch	=	26,
	gap	=	0
)
#dev.off()

t.test(DS.rates$DS3.D[1:9],DS.rates$DS5.D[1:9])