### USER INPUT: ###
directory <- "\\\\uniwa.uwa.edu.au\\userhome\\staff7\\00101127\\My Documents\\LLPS results\\20210617_Confocal\\plateII_B11_60x_0p10_2048.nd2_analysis(1-1000_th15_f2_1.2)"
no_of_drops <- 1029

# Reduce data (for each individual droplet, average out y values on equivalent
# x positions)?
  y_data_reduction <- TRUE

# Circle fit error (normalised RMSE) threshold for inclusion of droplet in
# final fit to find contact angle (relative error... e.g. 0.1 = 10%):
  RMSE_threshold <- 0.5


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
# f <- function(x) 100*(1/(sinh(x^0.7))) + 31


######  For fitting 3 coefficients: ###### 

f <- contact_angles ~ a*(1/(sinh(droplet_widths^b))) + c

model <- nls(formula = f, 
             data = results, 
             start = list(a=100, b=0.7, c=31))
summary(model)

a<-coef(model)[1]
b<-coef(model)[2]
c<-coef(model)[3]

func <- function(droplet_widths){ a*(1/(sinh(droplet_widths^b))) + c}

plt <- ggplot(data = results, mapping = aes(x=droplet_widths, y=contact_angles)) +
  geom_point(mapping = aes(colour=RMSE_norm)) +
  # geom_smooth() +
  stat_function(fun = func)
  # (+) coord_cartesian(ylim = c(25,100))

plt

# plotly::ggplotly(plt)

contact_angle_fit <- summary(model)$parameters[3,1]
contact_angle_fit_stderr <- summary(model)$parameters[3,2]


# 
# ###### For fitting only 2 coefficients: ######
# 
# f <- contact_angles ~ 100*(1/(sinh(droplet_widths^b))) + c
# 
# model <- nls(formula = f,
#              data = results,
#              start = list(b=0.7, c=31))
# 
# summary(model)
# 
# b<-coef(model)[1]
# c<-coef(model)[2]
# 
# func <- function(droplet_widths){ 100*(1/(sinh(droplet_widths^b))) + c}
# 
# plt <- ggplot(data = results, mapping = aes(droplet_widths, contact_angles)) +
#   geom_point() +
#   # geom_smooth() +
#   stat_function(fun = func)
# # (+) coord_cartesian(ylim = c(25,100))
# 
# plt
# 
# # plotly::ggplotly(plt)
# 
# contact_angle_fit <- summary(model)$parameters[2,1]
# contact_angle_fit_stderr <- summary(model)$parameters[2,2]
# 
# 
# 
# ###### For fitting only contact angle (1 coefficient): ######
# 
# a <- 100
# b <- 1.0
# 
# f <- contact_angles ~ a*(1/(sinh(droplet_widths^b))) + c
# 
# model <- nls(formula = f,
#              data = results,
#              start = list(c=31))
# 
# summary(model)
# 
# c<-coef(model)[1]
# 
# func <- function(droplet_widths){ a*(1/(sinh(droplet_widths^b))) + c}
# 
# plt <- ggplot(data = results, mapping = aes(droplet_widths, contact_angles)) +
#   geom_point(mapping = aes(colour=RMSE_norm)) +
#   # geom_smooth() +
#   stat_function(fun = func) +
#   ylim(0, 180)
# # (+) coord_cartesian(ylim = c(25,100))
# 
# plt
# 
# # plotly::ggplotly(plt)
# 
# contact_angle_fit <- summary(model)$parameters[1,1]
# contact_angle_fit_stderr <- summary(model)$parameters[1,2]


# ###### For fitting horizontal line ######


