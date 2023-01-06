##############################################################################################
## Activity curves for African Wild Dogs
##############################################################################################

# Load the packages

library("ggplot2")
library("dplyr")
library("glmmTMB")
library("splines2")


# Set the working directory // Update this to your working directory
setwd("C:/Users/kasim/Documents/1. My Documents/1. University/Activity and temperature/Proceedings of the Royal Society B/Reproducible data and code")

# Load pre-written functions
source("Scripts/Functions/predictPlot.R")

model_df_new2 <- read.csv('Data/activity_curves_dog.csv')

# Format needed columns
model_df_new2$status2 <- as.factor(model_df_new2$status2)
model_df_new2$hour_fc <- numFactor(model_df_new2$hour_fc)

#################################
## Regression modelling  
#################################

# Use periodic mSPlines
m_lux <- glmmTMB(status2 ~ mSpline(sunTime, df = 4, periodic = TRUE)*max_temp24Sc + season +
                   ou(hour_fc-1|date_fc) + (1|ID), 
                 family=binomial, data=model_df_new2, 
                 na.action = 'na.fail')

# Generate models
m_lux_aic <- MuMIn::dredge(m_lux, REML = FALSE)

