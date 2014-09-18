basefile	<-	'~/Dropbox/Shared Vito/cFP paper/Figures/'
#Code to automatically add PC compatibility
basefile2	<-	gsub('\\\\','/',path.expand(basefile))
xx			<-	try(setwd(basefile2))
if(class(xx)=="try-error"){basefile2 <- paste('C:/Users/Peter',substring(basefile,2),sep="")}
read.dir 	<- paste(basefile2,"Data used for figures (pulled by R)",sep="")
write.dir	<- paste(basefile2,"Figure parts",sep="")
setwd(read.dir)


dcga       <- read.csv('2014-04-10 erl BF DIP rate.csv')

# Get rid of DMSO past day 3
dcga    <- subset(dcga, !(Condition == 'DMSO' & Time_Day > 3) )

# Uniquely ID each set of data 
dcga$id <- paste(dcga$Expt,dcga$Well,dcga$Colony,dcga$Condition,dcga$Cell.line, sep="_")    
dcga$l2	<-	log2(dcga$Cell.count)
dcga[is.infinite(dcga$l2),'l2']	<-	0

# Log2 normalize to day 0
group	<- as.character(unique(dcga$id))
nl2		<-	vector()

for(i in group)	{
        temp	<- dcga[dcga$id == i,'l2'] -
        					dcga[dcga$id == i & dcga$Time_day == 0,'l2']
        ifelse(	match(i,group)==1, 
        		nl2 <-	temp,
        		nl2 <-	append(	nl2, temp))
        }
        
dcga$l2		<-	nl2        

rm(list=c('nl2','group','temp'))

source('estimate-slopes.R')