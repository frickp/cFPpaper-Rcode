# Load in cFP data for a variety of treatment conditions
cfp.data	<-	read.csv(textConnection(getURL(paste0(mybaseURL,'cFP-SingleAgentConcDep.csv'))))

### Code that scans the data frame, finds time points where all the cells disappear, 17Ad replace it
### with .999, so that the data can be displayed in log scale, yet easily found. 

cfp.data$Cell.count[cfp.data$Cell.count==0]	<-	0.999
cfp.data$log2 <- log2(cfp.data$Cell.count)


#List of all drugs to be combined with erlotinib
drugs	<-	c('CHX','FSK','TRM','X17A','SB','An','ABT','PLX')

#Add either DMSO 'D' or erlotinib 'E' to ID
cfp.data$ID1	<-	ifelse(cfp.data$Erl!=0,'E','D')

#Keep only non-zero conditions in ID
cfp.data$ID2	<-	apply(cfp.data[,drugs],MARGIN=c(1),FUN=function(x) 
{
	y	<-	subset(x,x>0)
	paste(names(y),y,sep='')
}
)

#Define unique IDs
cfp.data$ID3	<-	paste(cfp.data$ID1,cfp.data$ID2,sep='_')
#Define single-agent treatment
cfp.data[grep('character',cfp.data$ID3),'ID3']	<-	paste(substr(cfp.data[grep('character',
	cfp.data$ID3),'ID3'],1,2),'only',sep='')
#Final ID names
cfp.data$ID	<-	paste(
					paste(cfp.data$Expt,substr(cfp.data$Plate,1,2),sep='.'),
					paste(cfp.data$Well,cfp.data$Colony,sep=''),
					cfp.data$ID3,
					cfp.data$Cell.line,
				sep='_')
#Remove intermediate ID columns from data frame
cfp.data	<-	subset(cfp.data,select=names(cfp.data)[!names(cfp.data)%in%c('ID1','ID2','ID3')])

ID.list	<-	as.character(unique(cfp.data$ID))



###Normalize data to start at 0 from a log-scale
print('Processing and normalizing raw data...')
cfp.data$norm	<-	rep(0)
norm		<-	numeric()
for (i in ID.list)
	norm	<-	append(norm, 	subset(cfp.data, ID == i)$log2 - subset(cfp.data, ID == i)$log2[1])
	cfp.data$norm	<-	norm
print('Data normalization complete')

# Find All Slopes
print('Estimating best fit linear models for cell lineages...')
cfp.rates <- aggregate(cfp.data$log2, by=list(cfp.data$ID),
			FUN=function(x)
            {
				ifelse(length(x) < 5,time <-1:4,time <- 4:length(x))	#DMSO 3 days
                coef(lm(x[time] ~ time))[2] / 24						#estimated slope						 
			}			 
		)
colnames(cfp.rates) <- c('ID','rates')
