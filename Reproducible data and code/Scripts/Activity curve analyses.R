##############################################################################################
## Activity curves for cheetahs
##############################################################################################

#################################
## Set your environment
#################################

# Load the packages
library("ggplot2")
library("dplyr")
library("glmmTMB")
library("splines2")
library('ggpubr')

####
# User input required
####

# Set the working directory // Update this to your working directory
setwd("")

# Set your species: 'cheetah', 'leopard', 'dog', or 'lion'
species <- 'cheetah'

model_df_new2 <- read.csv(paste0('Data/activity_curves_', species, '.csv'))


####################################
# Formatting the data
####################################

####
## Set columns in correct format
#####

model_df_new2$status2 <- as.factor(model_df_new2$status2)
model_df_new2$hour_fc <- numFactor(model_df_new2$hour_fc)

####
## Split into dry and wet season
####

data_dry <- model_df_new2[model_df_new2$season == "dry", ]

data_wet <- model_df_new2[model_df_new2$season == "wet", ]

#################################
## Regression modelling
#################################

####
# Dry season only
####

m_dry <- glmmTMB(status2 ~ mSpline(sunTime, df = 4, periodic = TRUE)*max_temp24Sc +
                   ou(hour_fc-1|date_fc) + (1|ID),
                 family=binomial, data=data_dry,
                 na.action = 'na.fail')

m_dry_aic <- MuMIn::dredge(m_dry, REML = FALSE)

####
# Wet season only
####

m_wet <- glmmTMB(status2 ~ mSpline(sunTime, df = 4, periodic = TRUE)*max_temp24Sc +
                   ou(hour_fc-1|date_fc) + (1|ID),
                 family=binomial, data=data_wet,
                 na.action = 'na.fail')


m_wet_aic <- MuMIn::dredge(m_wet, REML = FALSE)



