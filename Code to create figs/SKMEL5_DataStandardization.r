d			<-	read.csv('BP_SKMEL5.csv')

##########################################################################################
#Step 1 is to rename the wells appropriately
##########################################################################################
#First get the row name then get the column name

#Get the col name
d$Col	<-	substr(d$Pre.well,10,11)
#Get the row name
for (i in 1:nrow(d))
{
	x	<-	d[i,'Pre.well']
	if(grepl('R02',x))	y <- 'B'
	else if(grepl('R03',x)) y <- 'C'
	else if(grepl('R04',x)) y <- 'D'
	else if(grepl('R05',x)) y <- 'E'
	else if(grepl('R06',x)) y <- 'F'
	else if(grepl('R07',x)) y <- 'G'
	
	d$Row[i]	<-	y
}

d$Well	<-	paste0(d$Row,d$Col)

##########################################################################################
#Step 2 is to delete excess columns and save the new csv file
##########################################################################################

d	<-	subset(d,select=-c(Pre.well,Row,Col))

SaveFile <- function(FileName = paste0(read.dir,'SKMEL5_cfp_final.csv'))
{
	write.csv(d,file=FileName)
}

#SaveFile()

