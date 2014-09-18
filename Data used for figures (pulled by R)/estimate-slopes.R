
# Find All Slopes
slopes <- aggregate(dcga$l2, by=list(dcga$Expt,dcga$Well,dcga$Colony,dcga$Condition,dcga$Cell.line),
                    FUN=function(x)
                    {
                        if(length(x) < 5) { 			# DMSO is only 0-3 days
                            time <- 1:4
                        #} else if (length(x) < 8) { 	# T790M is only 6 days
                        #    time <- 4:7
						#} else if (length(x) < 9) { 	# H1650 is only 7 days
                        #    time <- 4:8
                        } else if (length(x) > 11) {	# pretreatment data (-3:-1)
                        	time <- 7:14
                        } else { 						# Erlotinib otherwise is 10 days
                            time <- 4:length(x)
                        }
                        coef(lm(x[time] ~ time))[2]
                    })

colnames(slopes)	<-	c("Expt", "Well", "Colony", "Condition", "Cell.line", "Slope")

slopes$Slope		<-	slopes$Slope#/24 * log(2)
rates				<-	slopes
colnames(rates)[6]	<-	'Rate'
DIP.rates.dcga		<-	subset(rates, Condition == 'Erlotinib')