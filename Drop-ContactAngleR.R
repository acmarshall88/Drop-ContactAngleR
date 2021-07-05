### USER INPUT: ###
directory <- "\\\\uniwa.uwa.edu.au\\userhome\\staff7\\00101127\\My Documents\\LLPS results\\20210617_Confocal\\plateII_B11_60x_0p10_2048.nd2_analysis(1-1000_th15_f2_1.2)"
no_of_drops <- 1029

# Reduce data (for each individual droplet, average out y values on equivalent
# x positions)?
  y_data_reduction <- TRUE

# Circle fit error (normalised RMSE) threshold for inclusion of droplet in
# final fit to find contact angle (relative error... e.g. 0.1 = 10%):
  RMSE_threshold <- 0.3


#######################
droplet_widths <- c()
contact_angles <- c()
radii <- c() 

# Create Progress Bar for 'for' loop (below):
progress = txtProgressBar(min = 1, max = no_of_drops, initial = 1)

for (i in 1:no_of_drops) {

  # Write RMS errors of circlefit() to file (instead of printing to console):
  sink(file = "circfit RMS errors.csv", append = FALSE)
  
  xy_data_file <- paste(i, "_XZ_SurfacePx.csv", sep = "")
  
  xy <- read.csv(paste(directory, xy_data_file, sep = "\\"))
  
  # reduce data by averaging out y values on equivalent x positions:
  if (y_data_reduction==TRUE) {
    
    xy_averaged <- data.frame(
      X_um = unique(xy$X_micron),
      Ymean_um = NA
      )
    
    for (j in 1:nrow(xy_averaged)) {
      equiv_x_pos <- xy$X_micron == xy_averaged$X_um[j]
      xy_averaged$Ymean_um[j] <- mean(xy$Y_micron[equiv_x_pos])
    }
    
    #Assign cols that contain appropriate x and y coords:
    xdata <- xy_averaged$X_um
    ydata <- xy_averaged$Ymean_um
  
  } else {
    
    #Assign cols that contain appropriate x and y coords:
    xdata <- xy[,5]
    ydata <- xy[,6]
  }
  
  
  #Fit circle function to data (xy):
  cfit <- pracma::circlefit(xdata, ydata)
  radius <- cfit[3]
  # (append...)
  radii <- c(radii, radius)
  
  #translate data so that circle centre is (0,0):
  xy_trans <- data.frame(x = xdata-cfit[1], y = ydata-cfit[2])
  cfit_trans <- c(0, 0, radius)  
  ymin <- 0-cfit[2]
  
  # Function for a circle with centre (h,k) and radius r is: 
  # (x-h)^2 + (y-k)^2 = r^2 ... (#pythagorus)
  # If (h,k) = (0,0), 
  #   then  x^2 + y^2 = r^2
  #         x = sqrt(r^2 - y^2) (# 2 solns: same magnitude, opposite signs (-/+))
  
  # Coordinates for the points at either and of the drop-surface interface are
  # given by:
  xy_left <- c(-sqrt(cfit[3]^2 - ymin^2), ymin)
  xy_right <- c(sqrt(cfit[3]^2 - ymin^2), ymin)
  
  # therefore chord length is...
  droplet_width_temp <- xy_right[1]-xy_left[1]
  # (append...)
  droplet_widths <- c(droplet_widths, droplet_width_temp)
  
  #################################
  # xy_left and xy_right are euclidian vectors of length r (as circle centre is (0,0))
  
  # The angle(s) between the chord and tangent(s) between/at these 2 points on
  # the circle (i.e the Droplet Contact Angle) is equal to half the arc measure 
  # (i.e. angle between the 2 euclidian vectors)...
  
  # (from https://stackoverflow.com/questions/1897704/angle-between-two-vectors-in-r)
  #angle <- function(x,y){
    #dot.prod <- x%*%y 
    #norm.x <- norm(x,type="2")
    #norm.y <- norm(y,type="2")
    #theta <- acos(dot.prod / (norm.x * norm.y))
    #as.numeric(theta)
  #}
  # gives angle in radians. Muliply by 180/pi to get degrees...
  
  #ContactAngle <- 0.5 * (180/pi) * angle(xy_left, xy_right)
  ##################################
  
  # Alternatively, use isosceles triangle formed by radii and chord to calculate
  # angle (lambda) between chord and radius (basic trigonometry)...
  lambda <- acos((droplet_width_temp/2)/radius)
  
  # ... multiply by 180/pi to get degrees...
  lambda <- lambda * 180/pi
  
  if (ymin < 0) {
    lambda <- -1*lambda
  }
  
  # ... the angle(s) between the chord and tangent(s) (i.e the Droplet Contact 
  # Angle) is equal to 90 minus lambda...
  contact_angle_temp <- 90 - lambda
  # (append...)
  contact_angles <- c(contact_angles, contact_angle_temp)
  
  # (end output to file)
  sink()
  
  #Update progress bar:
  setTxtProgressBar(progress, i)
  
}



### visualisation ### (ONLY SHOWS FINAL DROPLET)

# circle fit:
library(conicfit)
cmodel <- calculateCircle(0, 0, radius, steps = 50, sector = c(0,180))
cmodel_df <- data.frame(x = cmodel[,1], y = cmodel[,2])

# make df of points of triangle representing surface of slide/plate and 
# associated arc measure:
slide_surface_arc <- data.frame(x = c(0, xy_left[1], (xy_right[1]), 0), 
                            y = c(0, xy_left[2], (xy_right[2]), 0))

library(ggplot2)
ggplot(NULL, aes(x, y)) +
  geom_point(data = xy_trans) +
  geom_path(data = cmodel_df) +
  geom_path(data = slide_surface_arc)

  

#####################

# Remove "RMS error: " label from RMS error list:
cfit_RMSEs <- read.csv("circfit RMS errors.csv", header = FALSE)
cfit_RMSEs <- tidyr::separate(data = cfit_RMSEs,
                         col = 1,
                         into = c(NA, NA, "cfit_RMSE"),
                         sep = " ")
# Convert to numeric:
cfit_RMSEs <- sapply(cfit_RMSEs, as.numeric)

# Normalise RMSE to radius of fitted circle (i.e. will now be %error):
RMSE_norm <- c(cfit_RMSEs/radii)

# Assign results to single data frame:
results <- data.frame(radii, cfit_RMSEs, RMSE_norm, droplet_widths, contact_angles)

# Remove rows with RMSE_norm greater than threshold value:
rows_to_keep <- RMSE_norm <= RMSE_threshold
results <- results[rows_to_keep,]


###### Fit hyperbolic function with 2 asymptotes: ######
# (vertical asymp = 0; horizontal asymp = contact angle)

# Manual approximation...
# f <- function(x) 80*(1/(sinh(x^0.6))) + 35
# (i.e. a=100, b=0.7, c=40)

# 'c' is the horizontal asymptote, and represents the theoretical contact angle
# of an infinitely large droplet.

# For fitting all 3 coefficients: # 

model_hyprblc <- nls(
              formula = contact_angles ~ a*(1/(sinh(droplet_widths^b))) + c, 
              data = results, 
              start = list(a=80, b=0.6, c=35),
              control = nls.control(maxiter = 200, minFactor = 1/4096))
summary(model_hyprblc)

a <- coef(model_hyprblc)[1]
b <- coef(model_hyprblc)[2]
c <- coef(model_hyprblc)[3]

func_hyprblc <- function(droplet_widths){ a*(1/(sinh(droplet_widths^b))) + c}

plt1 <- ggplot(data = results, mapping = aes(x=droplet_widths, y=contact_angles)) +
  geom_point(mapping = aes(colour=RMSE_norm)) +
  # geom_smooth() +
  stat_function(fun = func_hyprblc) +
  ylim(0,180)

plt1
# plotly::ggplotly(plt)

contact_angle_hyprblc <- summary(model_hyprblc)$parameters[3,1]
contact_angle_hyprblc_stderr <- summary(model_hyprblc)$parameters[3,2]


###### For fitting horizontal line ######
# (used when using a very strict RMSE_threshold... i.e. only large, well-fitted
# droplets should be included in the data)

model_linear <- nls(
              formula = contact_angles ~ y_int,
              data = results,
              start = list(y_int=35))
summary(model_linear)

y_int <- coef(model_linear)[1]

func_linear <- function(droplet_widths){y_int}

plt2 <- ggplot(data = results, mapping = aes(x=droplet_widths, y=contact_angles)) +
  geom_point(mapping = aes(colour=RMSE_norm)) +
  # geom_smooth() +
  stat_function(fun = func_linear) +
  ylim(0,180)

plt2
# plotly::ggplotly(plt)

contact_angle_linear <- summary(model_linear)$parameters[1,1]
contact_angle_linear_stderr <- summary(model_linear)$parameters[1,2]


###### For binning data before fitting hyperbolic function: ######

# set up dividers for bins:
breaks <- seq(0, 13, by = 0.5)

# specify interval/bin labels:
bin_num <- seq(0.5, 13, by = 0.5)
tags <- as.character(bin_num)

# assign values to bins:
bin <- cut(results$droplet_widths, 
           breaks=breaks, 
           include.lowest=TRUE, 
           right=FALSE, 
           labels=tags)
# inspect bins
summary(bin)

# Append bin column to results df:
results$bin <- bin

sorted_results <- results[order(results$bin),]

bin_mean_contact_angle <- c()

for (k in bin_num) {
  results_bin_logic <- sorted_results$bin == k
  results_single_bin <- na.omit(sorted_results[results_bin_logic,])
  bin_mean_contact_angle <- c(bin_mean_contact_angle, mean(results_single_bin$contact_angles))
}

binned_data <- data.frame(bin_num, bin_mean_contact_angle)

# Fit hyperbolic function: 

model_hyprblc_bins <- nls(formula = bin_mean_contact_angle ~ a*(1/(sinh(bin_num^b))) + c, 
             data = binned_data, 
             start = list(a=80, b=0.6, c=35),
             control = nls.control(maxiter = 200, minFactor = 1/4096))
summary(model_hyprblc_bins)

a<-coef(model_hyprblc_bins)[1]
b<-coef(model_hyprblc_bins)[2]
c<-coef(model_hyprblc_bins)[3]

# Visualisation #

func_hyprblc_bins <- function(bin_num){ a*(1/(sinh(bin_num^b))) + c}

plt3 <- ggplot() +
  geom_point(mapping=aes(x=droplet_widths, y=contact_angles, colour=RMSE_norm), data=results) +
  geom_point(mapping=aes(x=bin_num, y=bin_mean_contact_angle), data=binned_data, colour="magenta", size=5) +
  stat_function(fun=func_hyprblc_bins, colour="magenta", size=1) +
  ylim(0,180)

plt3
# plotly::ggplotly(plt)

contact_angle_hyprblc_bins <- summary(model_hyprblc_bins)$parameters[3,1]
contact_angle_hyprblc_bins_stderr <- summary(model_hyprblc_bins)$parameters[3,2]