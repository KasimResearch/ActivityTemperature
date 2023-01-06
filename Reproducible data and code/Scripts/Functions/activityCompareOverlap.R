##################################################
## QUANTIFY OVERLAP ACROSS THE TEMPERATURE RANGES
##################################################

###########
## SET ENVIRONMENT
###########

library(lubridate)
library("overlap")
library('dplyr')



###########
## WRITE OUT KEY FUNCTIONS
###########

# Generate random number of points (N), between a start (st), and end (et) time
timeDistribute <- function(N, st) {
       #st <- as.POSIXct(as.Date(st))
        st <- as.POSIXct(st, format = "%d/%m/%Y %H:%M:%S")
       et <- st + ((60*59)+59)
       dt <- as.numeric(difftime(et,st,unit="sec"))
       ev <- sort(runif(N, 0, dt))
       rt <- st + ev
       
       return(rt)
}

# Helper function to convert to radians
convert2rads <- function (distributedTimes) {
  
  
hours <- as.numeric(hour(distributedTimes))
minutes <- as.numeric(minute(distributedTimes))
seconds <- as.numeric(second(distributedTimes))

minHours <- minutes/60
secHours <- seconds/(60*60)

hourFrac <- hours + minHours + secHours

timeRad <- hourFrac * 0.261799

return(timeRad)

}

###########
## LOAD DATA 
###########

# # CHEETAH
# chMoveFull <- read.csv('Analyses/modelling_tables/hour/cheetah_movementFull_predicted.csv')
# 
# # Lion
# liDryMoveFull <- read.csv('Analyses/modelling_tables/hour/Lion_movementFull_dry_predicted.csv')
# 
# # CHEETAH
# awDenMoveFull <- read.csv('Analyses/modelling_tables/hour/dog_den_movement_full_predicted.csv')
# 
# # Lion
# awNodenMoveFull <- read.csv('Analyses/modelling_tables/hour/dog_noden_movement_full_predicted.csv')


## FOR TESTING

# speciesA <- awNodenMoveFull
# speciesB <- liDryMoveFull


activityCompareOverlap <- function(speciesA, speciesB) {

###########
## PREPARE FOR CREATING CURVES
###########

####
## Species A  
####  

speciesA <- speciesA %>%
  group_by(hour, max_temp24) %>%
             summarize(predictions = mean(fitted2))

speciesA$predictions <- speciesA$predictions*100


datalist <- list()

for (i in 1:nrow(speciesA)) {
  
  #Subset to row
  row <- speciesA[i, ]
  
  # Get the number of 'occurrences'
  nTimes <- as.numeric(round(row[1, 'predictions']))
  
  # Get the start time/hour
  st <- paste("01/01/1999 ", speciesA[i, 'hour'], ':00:00', sep = "")
  
  # Get your randomly distributed times
  times <- timeDistribute(nTimes, st)
  
  # Replicate row by the number of random distributions
  dataExpanded <- row[rep(seq_len(nrow(row)), nTimes), ]
  
  # Add the randomly distributed times
  dataExpanded$distributedTimes <- times
  
  datalist[[i]] <- dataExpanded
  
}

speciesA <- do.call(rbind, datalist)

# Convert to radians
speciesA$timeRad <- convert2rads(speciesA$distributedTimes)

speciesALow <- as.data.frame(speciesA[speciesA$max_temp24 == 20, ])
speciesAMedium <- as.data.frame(speciesA[speciesA$max_temp24 == 30, ])
speciesAHigh <- as.data.frame(speciesA[speciesA$max_temp24 == 40, ])

####
## Species B
####

speciesB <- speciesB %>%
  group_by(hour, max_temp24) %>%
  summarize(predictions = mean(fitted2))

speciesB$predictions <- speciesB$predictions*100


datalist <- list()

for (i in 1:nrow(speciesB)) {
  
  #Subset to row
  row <- speciesB[i, ]
  
  # Get the number of 'occurrences'
  nTimes <- as.numeric(round(row[1, 'predictions']))
  
  # Get the start time/hour
  st <- paste("01/01/1999 ", speciesB[i, 'hour'], ':00:00', sep = "")
  
  # Get your randomly distributed times
  times <- timeDistribute(nTimes, st)
  
  # Replicate row by the number of random distributions
  dataExpanded <- row[rep(seq_len(nrow(row)), nTimes), ]
  
  # Add the randomly distributed times
  dataExpanded$distributedTimes <- times
  
  datalist[[i]] <- dataExpanded
  
}

speciesB <- do.call(rbind, datalist)

# Convert to radians
speciesB$timeRad <- convert2rads(speciesB$distributedTimes)

speciesBLow <- speciesB[speciesB$max_temp24 == 20, ]
speciesBMedium <- speciesB[speciesB$max_temp24 == 30, ]
speciesBHigh <- speciesB[speciesB$max_temp24 == 40, ]


###########
## CREATE THE ACTIVITY CURVES
###########

#####
# Visualise the plots
#####
# 
# densityPlot(speciesALow$timeRad, rug=TRUE)
# densityPlot(speciesAMedium$timeRad, rug=TRUE)
# densityPlot(speciesAHigh$timeRad, rug=TRUE)
# 
# densityPlot(speciesBLow$timeRad, rug=TRUE)
# densityPlot(speciesBMedium$timeRad, rug=TRUE)
# densityPlot(speciesBHigh$timeRad, rug=TRUE)


####
## Bootstrap your overlap estimates
####

###
## LOW
###

# Calculate estimates of overlap:
overlapEstLow <- overlapEst(speciesALow$timeRad, speciesBLow$timeRad, 1000, type="Dhat4")

# Get boostraped values
bsLow <- bootstrap(speciesALow$timeRad, speciesBLow$timeRad, 1000, type="Dhat4") # takes a few seconds

# Get mean bootstraped value
bsMeanLow <- mean(overlapEstLow) 

# Get confidence intervals from bootstrapped values
ciLow <- bootCI(overlapEstLow, bsLow)['norm0', ]

summaryLow <- data.frame(mean = bsMeanLow, ciLow = ciLow[1], ciHigh = ciLow[2])
summaryLow$category <- 'low'
summaryLow$type <- 'overlapPackage'


###
## Medium
###

overlapEstMedium <- overlapEst(speciesAMedium$timeRad, speciesBMedium$timeRad, 10000, type="Dhat4")
bsMedium <- bootstrap(speciesAMedium$timeRad, speciesBMedium$timeRad, 1000, type="Dhat4") # takes a few seconds
bsMeanMedium <- mean(overlapEstMedium) 
ciMedium <- bootCI(overlapEstMedium, bsMedium)['norm0', ]

summaryMedium <- data.frame(mean = bsMeanMedium, ciLow = ciMedium[1], ciLow = ciMedium[2])
summaryMedium$category <- 'Medium'
summaryMedium$type <- 'overlapPackage'


###
## High
###

overlapEstHigh <- overlapEst(speciesAHigh$timeRad, speciesBHigh$timeRad, 1000, type="Dhat4")
bsHigh <- bootstrap(speciesAHigh$timeRad, speciesBHigh$timeRad, 1000, type="Dhat4") # takes a few seconds
bsMeanHigh <- mean(overlapEstHigh) 
ciHigh <- bootCI(overlapEstHigh, bsHigh)['norm0', ]

summaryHigh <- data.frame(mean = bsMeanHigh, ciLow = ciHigh[1], ciLow = ciHigh[2])
summaryHigh$category <- 'High'
summaryHigh$type <- 'overlapPackage'

###
## Prepare and bind together 
###

colnames(summaryLow) <- 
  colnames(summaryMedium) <- 
  colnames(summaryHigh) <- c('mean', 'ciLower', 'ciUpper', 'category', 'type')

output <- rbind(summaryLow, summaryMedium, summaryHigh)

rownames(output) <- NULL

return(output)

}



