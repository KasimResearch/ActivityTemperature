##################################################
## QUANTIFY OVERLAP ACROSS THE TEMPERATURE RANGES
##################################################

###########
## Set-up preamble
###########

# Set up environment
source('Scripts/Activity_overlaps_analyses_setup.R')

####
# User input required
####

# Set the working directory // Update this to your working directory
setwd(" ")

# Load cheetah data (change dry/wet as wanted)
cheetahData <- read.csv('Data/cheetah_activity_predictions_dry.csv')

# Load lion data (change dry/wet as wanted)
lionData <- read.csv('Data/lion_activity_predictions_dry.csv')

# Load AWD data (change dry/wet as wanted)
dogData <- read.csv('Data/dog_activity_predictions_dry.csv')

# Load leopard data (change dry/wet as wanted)
leopardData <- read.csv('Data/leopard_activity_predictions_dry.csv')

##################
## RUN FUNCTIONS
##################

#####
## Cheetah - lion
#####

# Create dataframes of percentile values to loop through (See above note on differences)
temp_percentile1 <- data.frame(percentile = 0:100, temperature = sort(unique(cheetahData$max_temp24Sc)))
temp_percentile2 <- data.frame(percentile = 0:100, temperature = sort(unique(lionData$max_temp24Sc)))

# Subset to evert fifth percentile for time constraints
temp_percentile1 <-  temp_percentile1[seq(1, nrow(temp_percentile1), by = 5), ]
temp_percentile2 <-  temp_percentile2[seq(1, nrow(temp_percentile2), by = 5), ]

# Create empty list to loop through results
datalist2 <- list()

startTime <- Sys.time()

# Loop through the dataframe
for (i in 1:nrow(temp_percentile1)) {

  loop_temp1 <- temp_percentile1[i, "temperature"]
  loop_temp2 <- temp_percentile2[i, "temperature"]
  loop_temp_percentile <- temp_percentile1[i, "percentile"]

  print(loop_temp_percentile)

  datalist2[[i]] <- activityCompareOverlapScaled(cheetahData, lionData, loop_temp1, loop_temp2, loop_temp_percentile)

  print(nrow(temp_percentile1)-i)

}

endTime <- Sys.time()

cheetah_lion <- do.call(rbind, datalist2)

rm(datalist2)


############
## Cheetah leopard
#############

# Create dataframes of percentile values to loop through (See above note on differences)
temp_percentile1 <- data.frame(percentile = 0:100, temperature = sort(unique(cheetahData$max_temp24Sc)))
temp_percentile2 <- data.frame(percentile = 0:100, temperature = sort(unique(leopardData$max_temp24Sc)))

# Subset to evert fifth percentile for time constraints
temp_percentile1 <-  temp_percentile1[seq(1, nrow(temp_percentile1), by = 5), ]
temp_percentile2 <-  temp_percentile2[seq(1, nrow(temp_percentile2), by = 5), ]

# Create empty list to loop through results
datalist2 <- list()

startTime <- Sys.time()

# Loop through the dataframe
for (i in 1:nrow(temp_percentile1)) {

  loop_temp1 <- temp_percentile1[i, "temperature"]
  loop_temp2 <- temp_percentile2[i, "temperature"]
  loop_temp_percentile <- temp_percentile1[i, "percentile"]

  print(loop_temp_percentile)

  datalist2[[i]] <- activityCompareOverlapScaled(cheetahData, leopardData, loop_temp1, loop_temp2, loop_temp_percentile)

  print(nrow(temp_percentile1)-i)

}

endTime <- Sys.time()

cheetah_leopard <- do.call(rbind, datalist2)

rm(datalist2)


############
## Cheetah dog
#############

# Create dataframes of percentile values to loop through (See above note on differences)
temp_percentile1 <- data.frame(percentile = 0:100, temperature = sort(unique(cheetahData$max_temp24Sc)))
temp_percentile2 <- data.frame(percentile = 0:100, temperature = sort(unique(dogData$max_temp24Sc)))

# Subset to evert fifth percentile for time constraints
temp_percentile1 <-  temp_percentile1[seq(1, nrow(temp_percentile1), by = 5), ]
temp_percentile2 <-  temp_percentile2[seq(1, nrow(temp_percentile2), by = 5), ]

# Create empty list to loop through results
datalist2 <- list()

startTime <- Sys.time()

# Loop through the dataframe
for (i in 1:nrow(temp_percentile1)) {

  loop_temp1 <- temp_percentile1[i, "temperature"]
  loop_temp2 <- temp_percentile2[i, "temperature"]
  loop_temp_percentile <- temp_percentile1[i, "percentile"]

  print(loop_temp_percentile)

  datalist2[[i]] <- activityCompareOverlapScaled(cheetahData, dogData, loop_temp1, loop_temp2, loop_temp_percentile)

  print(nrow(temp_percentile1)-i)

}

endTime <- Sys.time()

cheetah_dog <- do.call(rbind, datalist2)

plot(cheetahData[cheetahData$max_temp24Sc == temp_percentile1[i, "temperature"], "fitted2"], type = 'line')

rm(datalist2)

############
## leopard lion
#############

# Create dataframes of percentile values to loop through (See above note on differences)
temp_percentile1 <- data.frame(percentile = 0:100, temperature = sort(unique(leopardData$max_temp24Sc)))
temp_percentile2 <- data.frame(percentile = 0:100, temperature = sort(unique(lionData$max_temp24Sc)))

# Subset to evert fifth percentile for time constraints
temp_percentile1 <-  temp_percentile1[seq(1, nrow(temp_percentile1), by = 5), ]
temp_percentile2 <-  temp_percentile2[seq(1, nrow(temp_percentile2), by = 5), ]

# Create empty list to loop through results
datalist2 <- list()

startTime <- Sys.time()

# Loop through the dataframe
for (i in 1:nrow(temp_percentile1)) {

  loop_temp1 <- temp_percentile1[i, "temperature"]
  loop_temp2 <- temp_percentile2[i, "temperature"]
  loop_temp_percentile <- temp_percentile1[i, "percentile"]

  datalist2[[i]] <- activityCompareOverlapScaled(leopardData, lionData, loop_temp1, loop_temp2, loop_temp_percentile)

  print(nrow(temp_percentile1)-i)

}

endTime <- Sys.time()

leopard_lion <- do.call(rbind, datalist2)

rm(datalist2)

############
## leopard dog
#############

# Create dataframes of percentile values to loop through (See above note on differences)
temp_percentile1 <- data.frame(percentile = 0:100, temperature = sort(unique(leopardData$max_temp24Sc)))
temp_percentile2 <- data.frame(percentile = 0:100, temperature = sort(unique(dogData$max_temp24Sc)))

# Subset to evert fifth percentile for time constraints
temp_percentile1 <-  temp_percentile1[seq(1, nrow(temp_percentile1), by = 5), ]
temp_percentile2 <-  temp_percentile2[seq(1, nrow(temp_percentile2), by = 5), ]

# Create empty list to loop through results
datalist2 <- list()

startTime <- Sys.time()

# Loop through the dataframe
for (i in 1:nrow(temp_percentile1)) {

  loop_temp1 <- temp_percentile1[i, "temperature"]
  loop_temp2 <- temp_percentile2[i, "temperature"]
  loop_temp_percentile <- temp_percentile1[i, "percentile"]

  datalist2[[i]] <- activityCompareOverlapScaled(leopardData, dogData, loop_temp1, loop_temp2, loop_temp_percentile)

  print(nrow(temp_percentile1)-i)

}

endTime <- Sys.time()

leopard_dog <- do.call(rbind, datalist2)

rm(datalist2)

############
## lion-dog
#############

# Create dataframes of percentile values to loop through (See above note on differences)
temp_percentile1 <- data.frame(percentile = 0:100, temperature = sort(unique(lionData$max_temp24Sc)))
temp_percentile2 <- data.frame(percentile = 0:100, temperature = sort(unique(dogData$max_temp24Sc)))

# Subset to evert fifth percentile for time constraints
temp_percentile1 <-  temp_percentile1[seq(1, nrow(temp_percentile1), by = 5), ]
temp_percentile2 <-  temp_percentile2[seq(1, nrow(temp_percentile2), by = 5), ]

# Create empty list to loop through results
datalist2 <- list()

startTime <- Sys.time()

# Loop through the dataframe
for (i in 1:nrow(temp_percentile1)) {

  loop_temp1 <- temp_percentile1[i, "temperature"]
  loop_temp2 <- temp_percentile2[i, "temperature"]
  loop_temp_percentile <- temp_percentile1[i, "percentile"]

  datalist2[[i]] <- activityCompareOverlapScaled(lionData, dogData, loop_temp1, loop_temp2, loop_temp_percentile)

  print(nrow(temp_percentile1)-i)

}

endTime <- Sys.time()

leopard_dog <- do.call(rbind, datalist2)

rm(datalist2)


