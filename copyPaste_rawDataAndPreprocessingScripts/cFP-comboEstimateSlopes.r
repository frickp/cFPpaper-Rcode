
# Find All Slopes
#slopes <- aggregate(cfp$l2, by=list(cfp$Expt,cfp$Well,cfp$Colony,cfp$Condition,cfp$Cell.line),
a <- aggregate(cfp$l2, by=list(cfp$Well,cfp$Colony,cfp$Erl,cfp$CHX,cfp$FSK,cfp$TRM,cfp$X17A,cfp$Plate),
            FUN=function(x)
               {
					if(length(x) < 5) { 			# DMSO is only 0-3 days
                            time <- 1:4
                    } else { 						# otherwise is 10 days
                            time <- 4:length(x)
                    }
                    cbind(coef(lm(x[time] ~ time))[2],summary(lm(x[time] ~ time))$r.squared)
                })

slopes	<-	cbind(a[,1:8],a[9][,1][,1],a[9][,1][,2])
				
				
#colnames(slopes)	<-	c("Expt", "Well", "Colony", "Condition", "Cell.line", "Slope")
colnames(slopes)	<-	c("Well", "Colony","Erl", "CHX", "FSK","TRM","X17A","Plate","Slope","R2")

slopes$Slope		<-	slopes$Slope
rates				<-	slopes
colnames(rates)[6]	<-	'Rate'
#DIP.rates.cfp		<-	subset(rates, Condition == 'Erlotinib')




slopes$id	<-	paste(slopes$Well,slopes$Colony,
					paste("E",slopes$Erl,sep=""),
					paste("CHX",slopes$CHX,sep=""),
					paste("TRM",slopes$TRM,sep=""),
					paste("X17A",slopes$X17A,sep=""),
					paste("Plate",slopes$Plate,sep=""),sep="_"
				)