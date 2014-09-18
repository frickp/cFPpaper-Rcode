basefile	<-	'~/Dropbox/Lab/cFP paper/Figures/'
#Code to automatically add PC compatibility
basefile2	<-	gsub('\\\\','/',path.expand(basefile))
xx			<-	try(setwd(basefile2))
if(class(xx)=="try-error"){basefile2 <- paste('C:/Users/Peter',substring(basefile,2),sep="")}
read.dir 	<- paste(basefile2,"Data used for figures (pulled by R)",sep="")
setwd(read.dir)

#cfp       <- read.csv('2013-10-15 cFP E + CHX, FSK, TRM, 17A.csv')
cfp       <- read.csv('2013-10-15 cFP E + CHX, FSK, TRM, 17A v2.1.csv')
# Uniquely ID each set of data 
cfp$id <- paste(cfp$Well,cfp$Colony,paste("E",cfp$Erl,sep=""),paste("CHX",cfp$CHX,sep=""),paste("TRM",cfp$TRM,sep=""),paste("X17A",cfp$X17A,sep=""),paste("Plate",cfp$Plate,sep=""), sep="_")    
cfp$l2	<-	log2(cfp$Cell.count)
cfp[is.infinite(cfp$l2),'l2']	<-	0

# Log2 normalize to day 0
group	<- as.character(unique(cfp$id))
nl2		<-	vector()

for(i in group)	{
        temp	<- cfp[cfp$id == i,'l2'] -
        					cfp[cfp$id == i & cfp$Time.day == 0,'l2']
        ifelse(	match(i,group)==1, 
        		nl2 <-	temp,
        		nl2 <-	append(	nl2, temp))
        }
        
cfp$l2		<-	nl2        
rm(list=c('nl2','group','temp'))

