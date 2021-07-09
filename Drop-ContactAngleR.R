# README:
  # filenames for files containing droplet surface coordinates need to be in the
  # following format:
  # "[i]_[plane][file_suffix]", where... 
      # 'i' is an interger (droplet number)
      # 'plane' is the profile slice plane (i.e. XZ or XY)
      # 'file_suffix' must end in '.csv'
  # (e.g. 1_XZ_SurfacePx.csv)
  

######  USER INPUT:  ######
  
# (use "\\" for "\" when listing path to data)...
  directory <- "\\\\uniwa.uwa.edu.au\\userhome\\staff7\\00101127\\My Documents\\LLPS results\\20201127_Confocal\\ACM20201127\\plateII_H10_60x_0p175steps002.nd2_analysis_2-1000"
  no_of_drops <- 982

# Coordinate data file suffix:
  file_suffix <- "_SurfacePx.csv"

# List droplet profile slice planes to include (i.e. XZ and/or YZ):
  planes <- c("XZ", "YZ")
  
# Interpolate data (i.e. convert pixel coordinates to micron values before 
# circle-fitting)??
# (Note: If the voxel dimensions are equal (i.e. Vox is cubic), this will make
# no difference to the contact angle calculation)
  interpolate <- TRUE

# Reduce data in y-direction (for each individual droplet, average out y values
# on equivalent x positions)?
  y_data_reduction <- TRUE

# Minimum no. of points required to include drop in analysis (coordinate data 
# containing less points than this will be skipped):
# (Note: this is applied AFTER y-data reduction if TRUE)
  pixel_threshold <- 3
  
# Circle fit error (normalised RMSE) threshold for inclusion of droplet in
# final fit to find contact angle (relative error... e.g. 0.1 = 10%):
  RMSE_threshold <- 0.3

# Number of bins for assigning values to before fitting model (hyperbolic):
# (should be ~ 2-10% of no_of_drops)
  nbins <- no_of_drops*0.04

######  (END USER INPUT)  ######

  
######  INITIAL SET UP ######

closeAllConnections()
  
# Read in voxel dimensions for data interpolation:
if (interpolate==TRUE){
  vox_dimensions <- read.csv(paste(directory,"\\voxel_dimensions.csv",sep = ""))
  vox_width <- vox_dimensions$Vx_width
  vox_height <- vox_dimensions$Vx_height
  vox_depth <- vox_dimensions$Vx_depth
}

# Create empty vectors for filling with 'for' loop:
droplet_widths <- c()
contact_angles <- c()
radii <- c() 
skipped_droplets <- c()

# Create Progress Bar for 'for' loop (below):
progress = txtProgressBar(min = 1, max = no_of_drops, initial = 1)

# Remove previous file containing list of RMSEs (if exists):
if (file.exists("temp circfit RMS errors.csv")==TRUE) {
  unlink("temp circfit RMS errors.csv")
}

######  (END INITIAL SET UP) ######


######  START DATA PROCESSING ######
# Loop through csv files containing droplet coordinates, fitting circle to each:
for (l in planes) {
  cat("\nFitting circles & calculating contact angles for", l, "droplet profiles...\n")
  
  for (i in 1:no_of_drops) {
  
    # Write RMS errors of circlefit() to file (instead of printing to console):
    sink(file = "temp circfit RMS errors.csv", append = TRUE)
    
    xy_data_file <- paste(i, "_", l, file_suffix, sep = "")
    
    xy <- read.csv(paste(directory, xy_data_file, sep = "\\"))
    
    # Skip corrupted data:
    if (ncol(xy) < 3) {
      
      # (end output to file)
      sink()
      
      # Update progress bar:
      setTxtProgressBar(progress, i)
      
      # add droplet label to list
      skipped_droplets <- c(skipped_droplets, paste(i, "_", l, sep = ""))
      
      # Go back to top (next iteration)
      next
    }
    
    # Interpolate data (convert pixel coords to micron coords)
    if (interpolate==TRUE) {
      X_micron <- xy$X * vox_width
      Y_micron <- xy$Y * vox_depth
      xy <- cbind(xy, X_micron, Y_micron)
    }
    
    
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
      xdata <- xy$X_micron
      ydata <- xy$Y_micron
    }
    
    # Skip to next iteration if data are inappropriate for circle-fitting:
    if (length(xdata) < pixel_threshold | length(unique(ydata)) == 1) {
      
      # (end output to file)
      sink()
      
      # Update progress bar:
      setTxtProgressBar(progress, i)
      
      # add droplet label to list
      skipped_droplets <- c(skipped_droplets, paste(i, "_", l, sep = ""))
      
      # Go back to top (next iteration)
      next
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
}



### visualisation of circle-fitting (ONLY SHOWS FINAL DROPLET): ####

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

  

#### Prep data for model-fitting to find contact angle: ####

# Remove "RMS error: " label from RMS error list:
cfit_RMSEs <- read.csv("temp circfit RMS errors.csv", header = FALSE)
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

# Assign min and max values (for plots):
min_contact_angle <-min(results$contact_angles, na.rm = TRUE)
max_contact_angle <-max(results$contact_angles, na.rm = TRUE)
# min_droplet_width <-min(results$droplet_widths, na.rm = TRUE)
# max_droplet_width <-max(results$droplet_widths, na.rm = TRUE)
min_radii <-min(results$radii, na.rm = TRUE)
max_radii <-max(results$radii, na.rm = TRUE)



#### Fit horizontal line: ####
# (used when using a very strict data (pixel and RMSE) thresholds... 
# i.e. only large, well-fitted droplets should be included in the data)

model_linear <- nls(
  formula = contact_angles ~ y_int,
  data = results,
  start = list(y_int=35))
summary(model_linear)

contact_angle_linear <- summary(model_linear)$parameters[1,1]
contact_angle_linear_stderr <- summary(model_linear)$parameters[1,2]

y_int <- coef(model_linear)[1]

func_linear <- function(radii){y_int}

plt1 <- ggplot(data = results, mapping = aes(x=radii, y=contact_angles)) +
  labs(title = "horizontal line fit", 
       subtitle = paste("contact angle (y-intercept) =", 
                        round(contact_angle_linear),
                        "degrees",
                        sep = " "),
       x = "droplet radius (microns)",
       y = "contact angle (degrees)",
       colour = "RMS error\nof circle fit\n(relative)") +
    geom_point(mapping = aes(colour=RMSE_norm)) +
  # geom_smooth() +
  stat_function(fun = func_linear, colour = "dark blue", size=1) +
  ylim(0.9*min_contact_angle, 1.1*max_contact_angle) +
  xlim(0, 1.1*max_radii)

plt1
# plotly::ggplotly(plt)



#### Fit hyperbolic function (orthogonal asymptotes): ####
# (vertical asymp = 0; horizontal asymp = contact angle)

# Manual approximation...
# f <- function(x) 80*(1/(sinh(x^0.6))) + 35
# (i.e. a=80, b=0.7, c=40)

# 'c' is the horizontal asymptote, and represents the theoretical contact angle
# of an infinitely large droplet.

# For fitting all 3 coefficients: # 

model_hyprblc <- nls(
              formula = contact_angles ~ a/sinh(radii-vox_width/2)^b + c,
              data = results, 
              start = list(a=20, b=1, c=40),
              control = nls.control(maxiter = 200, minFactor = 1/4096))
summary(model_hyprblc)

contact_angle_hyprblc <- summary(model_hyprblc)$parameters[3,1]
contact_angle_hyprblc_stderr <- summary(model_hyprblc)$parameters[3,2]

# Visualisation: #

a_h <- coef(model_hyprblc)[1]
b_h <- coef(model_hyprblc)[2]
c_h <- coef(model_hyprblc)[3]

func_hyprblc <- function(radii){ a_h/sinh(radii-vox_width/2)^b_h + c_h}

plt2 <- ggplot(data = results, mapping = aes(x=radii, y=contact_angles)) +
  labs(title = "hyperbolic fit", 
       subtitle = paste("contact angle (asymptote) =", 
                        round(contact_angle_hyprblc),
                        "degrees",
                        sep = " "),
       x = "droplet radius (microns)",
       y = "contact angle (degrees)",
       colour = "RMS error\nof circle fit\n(relative)") +
  geom_point(mapping = aes(colour=RMSE_norm)) +
  stat_function(fun = func_hyprblc, colour = "dark blue", size=1) +
  ylim(0.9*min_contact_angle, 1.1*max_contact_angle) +
  xlim(0, 1.1*max_radii)

plt2
# plotly::ggplotly(plt)



#### Fit power function: ####

model_pwr <- nls(
  formula = contact_angles ~ a/(radii-vox_width/2)^b + c,
  data = results, 
  start = list(a=40, b=1, c=30),
  control = nls.control(maxiter = 200, minFactor = 1/4096))
summary(model_pwr)

contact_angle_pwr <- summary(model_pwr)$parameters[3,1]
contact_angle_pwr_stderr <- summary(model_pwr)$parameters[3,2]

# Visualisation: #

a_p <- coef(model_pwr)[1]
b_p <- coef(model_pwr)[2]
c_p <- coef(model_pwr)[3]

func_pwr <- function(radii){ a_p/(radii-vox_width/2)^b_p + c_p}

plt3 <- ggplot(data = results, mapping = aes(x=radii, y=contact_angles)) +
  labs(title = "power fit", 
       subtitle = paste("contact angle (asymptote) =", 
                        round(contact_angle_pwr),
                        "degrees",
                        sep = " "),
       x = "droplet radius (microns)",
       y = "contact angle (degrees)",
       colour = "RMS error\nof circle fit\n(relative)") +
  geom_point(mapping = aes(colour=RMSE_norm)) +
  stat_function(fun = func_pwr, colour = "dark blue", size=1) +
  ylim(0.9*min_contact_angle, 1.1*max_contact_angle) +
  xlim(0, 1.1*max_radii)

plt3
# plotly::ggplotly(plt)



###### For binning data before fitting: ######

# set up dividers for bins:
breaks <- seq(floor(min_radii), 
              ceiling(max_radii), 
              length.out = nbins+1)

# specify interval/bin labels:
bin_num <- breaks[2:length(breaks)]
labels <- as.character(bin_num)

# assign values to bins:
bin <- cut(results$radii, 
           breaks=breaks, 
           include.lowest=TRUE, 
           right=FALSE, 
           labels=labels)
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



#### Fit hyperbolic function (BINNED DATA): ####

model_hyprblc_bins <- nls(
            formula = bin_mean_contact_angle ~ a/sinh(bin_num-vox_width/2)^b + c, 
            data = binned_data, 
            start = list(a=25, b=0.15, c=30),
            control = nls.control(maxiter = 2000, minFactor = 1/(2^20)))
summary(model_hyprblc_bins)

contact_angle_hyprblc_bins <- summary(model_hyprblc_bins)$parameters[3,1]
contact_angle_hyprblc_bins_stderr <- summary(model_hyprblc_bins)$parameters[3,2]

# Visualisation: #

a_hbin<-coef(model_hyprblc_bins)[1]
b_hbin<-coef(model_hyprblc_bins)[2]
c_hbin<-coef(model_hyprblc_bins)[3]

func_hyprblc_bins <- function(bin_num){ a_hbin/sinh(bin_num-vox_width/2)^b_hbin + c_hbin}

plt4 <- ggplot() +
  labs(title = "hyperbolic fit to binned data", 
       subtitle = paste("contact angle (asymptote) =", 
                        round(contact_angle_hyprblc_bins),
                        "degrees",
                        sep = " "),
       x = "droplet radius (microns)",
       y = "contact angle (degrees)",
       colour = "RMS error\nof circle fit\n(relative)") +
  geom_point(mapping=aes(x=radii, y=contact_angles, colour=RMSE_norm), data=results) +
  geom_point(mapping=aes(x=bin_num, y=bin_mean_contact_angle), data=binned_data, colour="magenta", size=5) +
  stat_function(fun=func_hyprblc_bins, colour="magenta", size=1) +
  ylim(0.9*min_contact_angle, 1.1*max_contact_angle) +
  xlim(0, 1.1*max_radii)

plt4
# plotly::ggplotly(plt3)



#### Fit power function (BINNED DATA): ####

model_pwr_bins <- nls(formula = bin_mean_contact_angle ~ a/(bin_num-vox_width/2)^b + c, 
                          data = binned_data, 
                          start = list(a=40, b=1, c=30),
                          control = nls.control(maxiter = 2000, minFactor = 1/(2^20)))
summary(model_pwr_bins)

contact_angle_pwr_bins <- summary(model_pwr_bins)$parameters[3,1]
contact_angle_pwr_bins_stderr <- summary(model_pwr_bins)$parameters[3,2]

# Visualisation #

a_pbin<-coef(model_pwr_bins)[1]
b_pbin<-coef(model_pwr_bins)[2]
c_pbin<-coef(model_pwr_bins)[3]

func_pwr_bins <- function(bin_num){ a_hbin/(bin_num-vox_width/2)^b_hbin + c_hbin}

plt5 <- ggplot() +
  labs(title = "power fit to binned data", 
       subtitle = paste("contact angle (asymptote) =", 
                        round(contact_angle_pwr_bins),
                        "degrees",
                        sep = " "),
       x = "droplet radius (microns)",
       y = "contact angle (degrees)",
       colour = "RMS error\nof circle fit\n(relative)") +
  geom_point(mapping=aes(x=radii, y=contact_angles, colour=RMSE_norm), data=results) +
  geom_point(mapping=aes(x=bin_num, y=bin_mean_contact_angle), data=binned_data, colour="magenta", size=5) +
  stat_function(fun=func_pwr_bins, colour="magenta", size=1) +
  ylim(0.9*min_contact_angle, 1.1*max_contact_angle) +
  xlim(0, 1.1*max_radii)

plt5
# plotly::ggplotly(plt3)


