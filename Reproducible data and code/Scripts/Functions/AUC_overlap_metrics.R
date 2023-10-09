############
## CALCULATE THE AREA UNDER TWO CURVES
############

#####
## This function provides metrics of the area under two curves
#####

## INPUT ######
## Y1 is the predicted values of curve 1
## Y2 is the predicted values of curve 2
###############

################
## OUTPUT ######
##
## The function outputs a dataframe with one row and seven columns
##
## aucArea = the absolute area size under the two curves
## totalArea = the absolute area size under the max value of the two vectors at each point of x
## y1Area = the absolute area size under curve y1
## y2Area = the absolute area size under curve y2
##
## overProp = overall proportion of overlap between the two curves
## y1Prop = proportion of curve y1 which overlaps with curve y2
## y2Prop = proportion of curve y2 which overlaps with curve y1
##
###############


#########
## FOR TESTING
#########

# x <- seq(-2,2,0.001)
# y1 <- dnorm(x, 0, 1)
# y2 <- dnorm(x, -1,1)
# 
# 
# y2 <- speciesAHigh$predictions
# y1 <- speciesBHigh$predictions
# x <- 1:24
# 
# plot(x,y1, type = "l", bty = "l", lwd = 2, las = 1, ylab = "y")
# lines(x,y2, col = "red", lwd = 2)

#########

AUC_overlap_metrics <- function(y1, y2) {

####
## PREAMBLE
####

## Define the function to calculate the area from xy coordinates using
## a contour integral
area<-function(X){
  X<-rbind(X,X[1,])
  x<-X[,1]
  y<-X[,2] 
  lx<-length(x)
  abs(sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2)
}

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

####################
## FOR DIAGNOSTICS
####################

# Plot total area encompassing both curves
par(mfrow = c(1,2))
plot(x,y1, bty = "l", type = "n", las = 1, main = "Total area")
polygon(x_area, total, col = "purple")
lines(x,y1, col = "black", lwd = 5)
lines(x,y2, col = "green", lwd = 5)

# plot area under both curves
plot(x,y1, bty = "l", type = "n", las = 1, main = "Area under curves")
polygon(x_area, auc, col = "blue")
lines(x,y1, col = "black", lwd = 5)
lines(x,y2, col = "green", lwd = 5)

# plot area under curve 1
plot(x,y1, bty = "l", type = "n", las = 1, main = "Area under curves")
polygon(x_area, y1AUC, col = "blue")
lines(x,y1, col = "black", lwd = 5)
lines(x,y2, col = "green", lwd = 5)

# plot area under curve 1
plot(x,y1, bty = "l", type = "n", las = 1, main = "Area under curves")
polygon(x_area, y2AUC, col = "blue")
lines(x,y1, col = "black", lwd = 5)
lines(x,y2, col = "green", lwd = 5)

#################################

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

output <- data.frame(aucArea = auc_area, totalArea = total_area, y1Area = y1_area,
                      y2Area = y2_area, overProp = overProp, y1Prop = y1Prop,
                      y2Prop = y2Prop)

}


