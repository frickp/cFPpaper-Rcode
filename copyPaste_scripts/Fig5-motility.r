setwd('C:/Documents and Settings/BD/Desktop/Git/frickpl/motility 2011-07-15/spreadsheets')

library(stats4)
library(arm)

# Load all shape data
wells <- c('C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09')
data <- list()
for(w in wells)
{
  cat("Reading ", w, '\n')
  data[[w]] <- read.csv(paste('2011-07-15_000Well',w,'_shapes.csv', sep=''))
}


# Set experimental conditions
sample_time       <- 12 # In minutes
samples           <- 20 
microns_per_pixel <- 2*0.312556

long_enough <- function(x, length)
{
  result <- NULL
  for(i in unique(x$Cell.ID))
  {
    d <- subset(x, Cell.ID==i)
    if(length(d$Centroid.1) >= length)  # Are there twenty observations or more for a track?
    { 
      result <- append(result, i)
    }
  }
  result
}

# Detail Results
single.cell.speed <- function(x, minutes, microns_per_pixel, length=20)
{
  result <- NULL
  ids    <- NULL
  for(i in unique(x$Cell.ID))
  {
    d <- subset(x, Cell.ID==i)
    if(length(d$Centroid.1) >= length)  # Are there twenty observations or more for a track?
    {
      coords <- complex(d$Centroid.1, d$Centroid.2)
      s <- Mod(coords[2:(length(coords))] - coords[1:(length(coords)-1)]) * microns_per_pixel
      result <- append(result, mean(s))
      ids <- append(ids, i)
    }
  }
  list(id = ids, speed = result/minutes * 60)
}
result <- NULL
for(w in wells)
{
  x <- single.cell.speed(data[[w]], sample_time, microns_per_pixel, samples)
  d <- data.frame(well=rep(w, length(x$id)), id=x$id, speed=x$speed)
  result <- rbind(result, d)
}

result$well <- factor(result$well, levels=c('C02', 'D02', 'C03', 'D03', 'C04', 'D04', 'C05', 'D05', 'C06', 'D06', 'C07', 'D07', 'C08', 'D08', 'C09', 'D09'), ordered=TRUE)
result$logspeed <- log(result$speed)
a <- aov(logspeed ~ well, result)
summary(a)                                    #show the summary table
print(model.tables(a,"means"),digits=3)

#report the means and the number of subjects/cell

DS3	<-	append(subset(result,well=="C04")$logspeed,subset(result,well=="D04")$logspeed)
DS5	<-	append(subset(result,well=="C06")$logspeed,subset(result,well=="D06")$logspeed)

dev.new(width=2,height=4)
bp	<-	paste(write.dir,"/Fig XX motility.pdf",sep="")
#pdf(bp,width=2, height=4)
par(font.lab=2,las=2)
boxplot(DS3,DS5,names=c("DS3","DS5"),notch=T,col='lightgrey',ylab="log(microns/hour)",ylim=c(0,5.5))
t.test(DS3,DS5)
#dev.off()