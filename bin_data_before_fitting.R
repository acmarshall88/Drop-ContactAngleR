# set up cut-off values 
breaks <- seq(0, 13, by = 0.5)

# specify interval/bin labels
tags <- seq(0.5, 13, by = 0.5)
tags <- as.character(tags)

# bucketing values into bins
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

bin_num <- seq(0.5, 13, by = 0.5)

bin_mean_contact_angle <- c()

for (k in bin_num) {
   results_bin_logic <- sorted_results$bin == k
   results_single_bin <- na.omit(sorted_results[results_bin_logic,])
   bin_mean_contact_angle <- c(bin_mean_contact_angle, mean(results_single_bin$contact_angles))
}

bin_mean_contact_angle

binned_data <- data.frame(bin_num, bin_mean_contact_angle)

#### Visualisation ####

f <- bin_mean_contact_angle ~ a*(1/(sinh(bin_num^b))) + c

model <- nls(formula = f, 
             data = binned_data, 
             start = list(a=80, b=0.6, c=38),
             control = nls.control(maxiter = 200, minFactor = 1/4096))
summary(model)

a<-coef(model)[1]
b<-coef(model)[2]
c<-coef(model)[3]

func <- function(bin_num){ a*(1/(sinh(bin_num^b))) + c}

# func <- function(droplet_widths){ 0.3*(1/(sinh(droplet_widths^1.2))) + 0.02}


plt <- ggplot(data = binned_data, mapping = aes(x=bin_num, y=bin_mean_contact_angle)) +
   geom_point() +
   # geom_point(mapping = aes(colour=RMSE_norm))
   # geom_smooth() +
   stat_function(fun = func)
   # ylim(0,0.12)
# (+) coord_cartesian(ylim = c(25,100))

plt

# plotly::ggplotly(plt)

