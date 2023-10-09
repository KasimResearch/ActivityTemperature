
####################################################################
## Activity overlap analyses setup
####################################################################

###########
## SET ENVIRONMENT
###########

library('lubridate')
library("overlap")
library('boot')
library('dplyr')
library('ggplot2')
library('gridExtra')
library('ggpubr')


#########################
# DEFINE KEY FUNCTIONS
#########################

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

convert2Rads <- function(timestamp) {

  hours <- as.numeric(hour(timestamp))
  minutes <- as.numeric(minute(timestamp))
  seconds <- as.numeric(second(timestamp))

  minHours <- minutes/60
  secHours <- seconds/(60*60)

  hourFrac <- hours + minHours + secHours

  timeRadians <- hourFrac * 0.261799

  return(timeRadians)

}

source('Scripts/Functions/bhatt.coef.R')


# Function to bootstrap (manually, used when data1 and data2 vectors differ in length)
bootBA <- function(data1, data2, bootLength = 1000) {

  datalist <- list()

  for (i in 1:bootLength) {

    # Resample data
    resample1 <- sample(data1, replace = T)
    resample2 <- sample(data2, replace = T)

    datalist[[i]] <- bhatt.coeff(resample1, resample2)


  }

  bootData <- (do.call(rbind, datalist))

  hist(bootData)

  return(data.frame(mean = mean(bootData),
                    ciLow = as.numeric(quantile(bootData, c(0.05,0.95)))[1],
                    ciHigh = as.numeric(quantile(bootData, c(0.05,0.95)))[2]))

}

## This function calculates the absolute area under two curves
## Y1 is the predicted values of curve 1
## Y2 is the predicted values of curve 2
AUC2 <- function (y1, y2) {

  x <- 1:length(y1)
  # vectors to store the area under the curves and the total
  # under = min of the two vlues
  # total = max of the two values
  auc <- rep(NA, length(x))
  total <- rep(NA, length(x))
  for(i in 1:length(auc)){
    auc[i] <- min(c(y1[i], y2[i]))
    total[i] <- max(c(y1[i], y2[i]))
  }
  # add zeroes to start and end of both vectors
  # so that the y values and x values meet at the
  # bottom
  auc <- c(0,auc,0)
  total <- c(0, total, 0)
  # repeat the min and max of x for this reason
  x_area <- c(min(x), x, max(x))

  # a function to calculate area from xy coords
  # via a contour integral.
  area<-function(X){
    X<-rbind(X,X[1,])
    x<-X[,1]
    y<-X[,2]
    lx<-length(x)
    abs(sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2)
  }

  # area under the two curves
  auc_area <- area(cbind(x_area, auc))

  return (auc_area)

}

## This function calculates the proportion of y1 under the y2 curve
AUC3 <- function (y1, y2) {

  x <- 1:length(y1)

  #########
  ## CREATING VECTORS TO STORE AUC VALUES
  #########

  # y1AUC = AUC for y1
  # y2AUC = AUC for y2
  # under = AUC for min of the two vectors at each point of x
  # total = AUC for max of the two vectors at each point of x

  # Create initial vectors
  y1AUC <- rep(NA, length(x))
  y2AUC <- rep(NA, length(x))
  auc <- rep(NA, length(x))
  total <- rep(NA, length(x))

  for(i in 1:length(auc)){
    y1AUC[i] <- y1[i] # redundant way to write it, but it's for clarity in my head
    y2AUC[i] <- y2[i] # redundant way to write it, but it's for clarity in my head
    auc[i] <- min(c(y1[i], y2[i]))
    total[i] <- max(c(y1[i], y2[i]))
  }

  # Add zeroes to start and end of the vectors so that the y values and x values
  # meet at the bottom
  y1AUC <- c(0,y1AUC,0)
  y2AUC <- c(0, y2AUC, 0)
  auc <- c(0,auc,0)
  total <- c(0, total, 0)

  # Repeat the min and max of x for the above reason
  x_area <- c(min(x), x, max(x))
  #################################

  # a function to calculate area from xy coords
  # via a contour integral.
  area<-function(X){
    X<-rbind(X,X[1,])
    x<-X[,1]
    y<-X[,2]
    lx<-length(x)
    abs(sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2)
  }

  ######
  ## CALCULATE THE AREA UNDER THE TWO CURVES
  ######

  # Area under the two curves
  auc_area <- area(cbind(x_area, auc))

  # Total area
  total_area <- area(cbind(x_area, total))

  # Area under curve 1
  y1_area <- area(cbind(x_area, y1AUC))

  # Area under curve 2
  y2_area <- area(cbind(x_area, y2AUC))

  #####
  ## WORK OUT PROPORTIONS
  #####

  ####
  # Prop. of Y1 overlapping
  ####

  # Overall proportion overlap
  overProp <- auc_area / total_area

  # Y1 proportion overlap
  y1Prop <-  auc_area / y1_area

  # Y1 proportion overlap
  y2Prop <-  auc_area / y2_area

  return(y1Prop)

}

## This function calculates the proportion of y1 under the y2 curve
AUC4 <- function (y1, y2) {

  x <- 1:length(y1)

  #########
  ## CREATING VECTORS TO STORE AUC VALUES
  #########

  # y1AUC = AUC for y1
  # y2AUC = AUC for y2
  # under = AUC for min of the two vectors at each point of x
  # total = AUC for max of the two vectors at each point of x

  # Create initial vectors
  y1AUC <- rep(NA, length(x))
  y2AUC <- rep(NA, length(x))
  auc <- rep(NA, length(x))
  total <- rep(NA, length(x))

  for(i in 1:length(auc)){
    y1AUC[i] <- y1[i] # redundant way to write it, but it's for clarity in my head
    y2AUC[i] <- y2[i] # redundant way to write it, but it's for clarity in my head
    auc[i] <- min(c(y1[i], y2[i]))
    total[i] <- max(c(y1[i], y2[i]))
  }

  # Add zeroes to start and end of the vectors so that the y values and x values
  # meet at the bottom
  y1AUC <- c(0,y1AUC,0)
  y2AUC <- c(0, y2AUC, 0)
  auc <- c(0,auc,0)
  total <- c(0, total, 0)

  # Repeat the min and max of x for the above reason
  x_area <- c(min(x), x, max(x))


  #################################

  # a function to calculate area from xy coords
  # via a contour integral.
  area<-function(X){
    X<-rbind(X,X[1,])
    x<-X[,1]
    y<-X[,2]
    lx<-length(x)
    abs(sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2)
  }

  ######
  ## CALCULATE THE AREA UNDER THE TWO CURVES
  ######

  # Area under the two curves
  auc_area <- area(cbind(x_area, auc))

  # Total area
  total_area <- area(cbind(x_area, total))

  # Area under curve 1
  y1_area <- area(cbind(x_area, y1AUC))

  # Area under curve 2
  y2_area <- area(cbind(x_area, y2AUC))

  #####
  ## WORK OUT PROPORTIONS
  #####

  ####
  # Prop. of Y1 overlapping
  ####

  # Overall proportion overlap
  overProp <- auc_area / total_area

  # Y1 proportion overlap
  y1Prop <-  auc_area / y1_area

  # Y1 proportion overlap
  y2Prop <-  auc_area / y2_area

  return(y2Prop)

}

## This function calculates the proportion of total overlap (auc area/total area)
AUC5 <- function (y1, y2) {

  x <- 1:length(y1)

  #########
  ## CREATING VECTORS TO STORE AUC VALUES
  #########

  # y1AUC = AUC for y1
  # y2AUC = AUC for y2
  # under = AUC for min of the two vectors at each point of x
  # total = AUC for max of the two vectors at each point of x

  # Create initial vectors
  y1AUC <- rep(NA, length(x))
  y2AUC <- rep(NA, length(x))
  auc <- rep(NA, length(x))
  total <- rep(NA, length(x))

  for(i in 1:length(auc)){
    y1AUC[i] <- y1[i] # redundant way to write it, but it's for clarity in my head
    y2AUC[i] <- y2[i] # redundant way to write it, but it's for clarity in my head
    auc[i] <- min(c(y1[i], y2[i]))
    total[i] <- max(c(y1[i], y2[i]))
  }

  # Add zeroes to start and end of the vectors so that the y values and x values
  # meet at the bottom
  y1AUC <- c(0,y1AUC,0)
  y2AUC <- c(0, y2AUC, 0)
  auc <- c(0,auc,0)
  total <- c(0, total, 0)

  # Repeat the min and max of x for the above reason
  x_area <- c(min(x), x, max(x))


  #################################

  # a function to calculate area from xy coords
  # via a contour integral.
  area<-function(X){
    X<-rbind(X,X[1,])
    x<-X[,1]
    y<-X[,2]
    lx<-length(x)
    abs(sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2)
  }

  ######
  ## CALCULATE THE AREA UNDER THE TWO CURVES
  ######

  # Area under the two curves
  auc_area <- area(cbind(x_area, auc))

  # Total area
  total_area <- area(cbind(x_area, total))

  # Area under curve 1
  y1_area <- area(cbind(x_area, y1AUC))

  # Area under curve 2
  y2_area <- area(cbind(x_area, y2AUC))

  #####
  ## WORK OUT PROPORTIONS
  #####

  ####
  # Prop. of Y1 overlapping
  ####

  # Overall proportion overlap
  overProp <- auc_area / total_area

  # Y1 proportion overlap
  y1Prop <-  auc_area / y1_area

  # Y1 proportion overlap
  y2Prop <-  auc_area / y2_area

  return(overProp)

}

# Takes in two equal vectors (x and y) and bootstraps them i times
# Returns BA coefficient and BCa confidence intervals for them
bootBootBA <- function (x, y, i = 1000) {# x = vector 1, y = vector 2, i = number of boot iterations

  # move data into data frame
  df <- data.frame(x = x, y = y)

  # define statistic to be bootstrapped
  # use boot's built-in `i` index to resample vector one
  # create a custom resample index to resample vector two
  BA2 <- function(df, i) {
    bhatt.coeff(df$x[i], df$y[sample(1:nrow(df), nrow(df), TRUE)])
  }

  # run bootstrap
  boot.output <- boot(df, statistic = BA2, R = i)

  # calculate CIs
  limits <- boot.ci(boot.output)

  output <- data.frame(mean = limits$t0, ciLow = limits$bca[4],
                       ciHigh = limits$bca[5])

  return(output)

}

# Takes in two equal vectors (x and y) and bootstraps them i times
# Returns AUC and confidence intervals for them
bootBootAUC <- function (x, y, i = 1000) {# x = vector 1, y = vector 2, i = number of boot iterations

  # move data into data frame
  df <- data.frame(x = x, y = y)

  # define statistic to be bootstrapped
  # use boot's built-in `i` index to resample vector one
  # create a custom resample index to resample vector two
  AUCboot <- function(df, i) {
    AUC2(df$x[i], df$y[sample(1:nrow(df), nrow(df), TRUE)])
  }

  # run bootstrap
  boot.output <- boot(df, statistic = AUCboot, R = i)

  # calculate CIs
  limits <- boot.ci(boot.output)

  output <- data.frame(mean = limits$t0, ciLow = limits$bca[4],
                       ciHigh = limits$bca[5])

  return(output)

}

bootBootAUC3 <- function (x, y, i = 1000) {# x = vector 1, y = vector 2, i = number of boot iterations

  # move data into data frame
  df <- data.frame(x = x, y = y)

  # define statistic to be bootstrapped
  # use boot's built-in `i` index to resample vector one
  # create a custom resample index to resample vector two
  AUCboot <- function(df, i) {
    AUC3(df$x[i], df$y[sample(1:nrow(df), nrow(df), TRUE)])
  }

  # run bootstrap
  boot.output <- boot(df, statistic = AUCboot, R = i)

  # calculate CIs
  limits <- boot.ci(boot.output)

  output <- data.frame(mean = limits$t0, ciLow = limits$bca[4],
                       ciHigh = limits$bca[5])

  return(output)

}

bootBootAUC4 <- function (x, y, i = 1000) {# x = vector 1, y = vector 2, i = number of boot iterations

  # move data into data frame
  df <- data.frame(x = x, y = y)

  # define statistic to be bootstrapped
  # use boot's built-in `i` index to resample vector one
  # create a custom resample index to resample vector two
  AUCboot <- function(df, i) {
    AUC4(df$x[i], df$y[sample(1:nrow(df), nrow(df), TRUE)])
  }

  # run bootstrap
  boot.output <- boot(df, statistic = AUCboot, R = i)

  # calculate CIs
  limits <- boot.ci(boot.output)

  output <- data.frame(mean = limits$t0, ciLow = limits$bca[4],
                       ciHigh = limits$bca[5])

  return(output)

}

bootBootAUC5 <- function (x, y, i = 1000) {# x = vector 1, y = vector 2, i = number of boot iterations

  # move data into data frame
  df <- data.frame(x = x, y = y)

  # define statistic to be bootstrapped
  # use boot's built-in `i` index to resample vector one
  # create a custom resample index to resample vector two
  AUCboot <- function(df, i) {
    AUC5(df$x[i], df$y[sample(1:nrow(df), nrow(df), TRUE)])
  }

  # run bootstrap
  boot.output <- boot(df, statistic = AUCboot, R = i)

  # calculate CIs
  limits <- boot.ci(boot.output)

  output <- data.frame(mean = limits$t0, ciLow = limits$bca[4],
                       ciHigh = limits$bca[5])

  return(output)

}

activityCompare <- function(speciesA1, speciesB1, iterations = 1000000) {

  ###########
  ## PREPARE FOR CREATING CURVES
  ###########

  ## Species A

  speciesA <- speciesA1 %>%
    group_by(sunTime, max_temp24Sc) %>% # Changed
    summarize(predictions = mean(fitted2))


  ## Species B
  speciesB <- speciesB1 %>%
    group_by(sunTime, max_temp24Sc) %>% # Changed
    summarize(predictions = mean(fitted2))



  ####
  # CALCULATE OVERLAP PACKAGE METRIC
  ####

  summaryOverlap <- activityCompareOverlapScaled(speciesA1, speciesB1)

  # colnames(summaryBA) <-
  #   colnames(summaryAUC) <-
  #   colnames(summaryAUC3) <-
  #   colnames(summaryAUC4) <-
  #   colnames(summaryAUC5) <-
  #   colnames(summaryOverlap)
  #
  #
  #
  # output <- rbind(summaryBA, summaryAUC, summaryAUC3, summaryAUC4, summaryAUC5, summaryOverlap) #
  #
  # return(output)

  return(summaryOverlap)

}


activityComparePlots <- function(speciesA, speciesB,
                                 speciesALabel = "speciesA",
                                 speciesBLabel = "speciesB") {


  ###########
  ## PREPARE FOR CREATING CURVES
  ###########

  ## Species A

  speciesA <- speciesA %>%
    group_by(sunTime, max_temp24Sc) %>%
    summarize(predictions = mean(fitted2))

  speciesALow <- speciesA[speciesA$max_temp24Sc == -2, ]
  speciesAHigh <- speciesA[speciesA$max_temp24Sc == 2, ]


  ## Species B
  speciesB <- speciesB %>%
    group_by(sunTime, max_temp24Sc) %>%
    summarize(predictions = mean(fitted2))

  speciesBLow <- speciesB[speciesB$max_temp24Sc == -2, ]
  speciesBHigh <- speciesB[speciesB$max_temp24Sc == 2, ]

  # PLOTS

  output <- list()

  # speciesALabel <- deparse(substitute(speciesA))
  # speciesBLabel <- deparse(substitute(speciesB))

  output[[1]] <- ggplot() +
    geom_smooth(data = speciesALow, aes(x = sunTime, y = predictions, linetype = 'species A'),
                span = 0.2, se = FALSE) +
    geom_smooth(data = speciesBLow, aes(x = sunTime, y = predictions, linetype = 'species B'),
                span = 0.2, se = FALSE) +
    ggtitle('Low temperature overlap') +
    scale_linetype_manual(values = c(1, 2),
                          labels = c(speciesALabel, speciesBLabel)) +
    scale_color_manual(values = c("black", "red"),
                       labels = c(speciesALabel, speciesBLabel)) +
    theme(legend.position="none") +
    scale_y_continuous(limits=c(0, 1)) +
    geom_vline(xintercept = 1.570796, linetype="dotted",
               color = "black") +
    geom_vline(xintercept = 4.712389, linetype="dotted",
               color = "black") +
    geom_text(aes(x = 1.9, label =  "sunrise", y = 0.75),
              colour = "black",
              text = element_text(size=10)) +
    geom_text(aes(x = 5.05, label =  "sunset", y = 0.75),
              colour = "black",
              text = element_text(size=10))




  output[[2]] <- ggplot() +
    geom_smooth(data = speciesAHigh, aes(x = sunTime, y = predictions, linetype = 'species A'),
                span = 0.2, se = FALSE) +
    geom_smooth(data = speciesBHigh, aes(x = sunTime, y = predictions, linetype = 'species B'),
                span = 0.2, se = FALSE) +
    ggtitle('High temperature overlap') +
    scale_linetype_manual(values = c(1, 2),
                          labels = c(speciesALabel, speciesBLabel)) +
    scale_color_manual(values = c("black", "red"),
                       labels = c(speciesALabel, speciesBLabel)) +
    theme(legend.position="none") +
    scale_y_continuous(limits=c(0, 1)) +
    geom_vline(xintercept = 1.570796, linetype="dotted",
               color = "black") +
    geom_vline(xintercept = 4.712389, linetype="dotted",
               color = "black") +
    geom_text(aes(x = 1.9, label =  "sunrise", y = 0.75),
              colour = "black",
              text = element_text(size=10)) +
    geom_text(aes(x = 5.05, label =  "sunset", y = 0.75),
              colour = "black",
              text = element_text(size=10))



  return(output)

}


#######
## LOAD FUNCTIONS
#######

source("Scripts/Functions/AUC_overlap_metrics.R")
source("Scripts/Functions/activityCompareOverlap.R")
source("Scripts/Functions/activityCompareOverlapScaled.R")
source("Scripts/Functions/activityCompareOverlapProps.R")
