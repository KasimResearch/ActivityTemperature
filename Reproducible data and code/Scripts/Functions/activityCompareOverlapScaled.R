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


activityCompareOverlapScaled <- function(speciesA, speciesB, temperature1, temperature2, temperature_percentile) {

  ###########
  ## PREPARE FOR CREATING CURVES
  ###########

  ####
  ## Species A
  ####

  speciesA <- speciesA %>%
    group_by(hour, max_temp24Sc) %>%
    summarize(predictions = mean(fitted2))

  speciesA$predictions <- speciesA$predictions*100

  # Subset to specific temperature
  speciesA <- as.data.frame(speciesA[speciesA$max_temp24Sc == temperature1, ])

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

  #speciesALow <- as.data.frame(speciesA[speciesA$max_temp24Sc == temperature, ])

  # Redundant renaming, but needing to keep formatting from legacy code
  speciesALow <- speciesA

  ####
  ## Species B
  ####

  speciesB <- speciesB %>%
    group_by(hour, max_temp24Sc) %>%
    summarize(predictions = mean(fitted2))

  speciesB$predictions <- speciesB$predictions*100

  # Subset to specific temperature
  speciesB <- as.data.frame(speciesB[speciesB$max_temp24Sc == temperature2, ])

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

  #speciesBLow <- speciesB[speciesB$max_temp24Sc == temperature, ]

  speciesBLow <- speciesB


  ###########
  ## CREATE THE ACTIVITY CURVES
  ###########

  #####
  # Visualise the plots
  #####
  #
#  densityPlot(speciesALow$timeRad, rug=TRUE)
 # densityPlot(speciesBLow$timeRad, rug=TRUE)


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
  summaryLow$category <- temperature_percentile
  summaryLow$type <- 'overlapPackage'

  ###
  ## Prepare and bind together
  ###

  colnames(summaryLow) <- c('mean', 'ciLower', 'ciUpper', 'temperature_percentile', 'type')

  output <- summaryLow

  rownames(output) <- NULL

  return(output)

}



