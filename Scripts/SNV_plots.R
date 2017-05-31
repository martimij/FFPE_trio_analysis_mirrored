# Martina Mijuskovic
# FFPE project
# SNV comparison and correlation plots
# May 2017


######### Load libraries ######### 

# NOTE: Installing R packages on HPC: use lib = "~/R/x86_64-pc-linux-gnu-library/3.X"

library(dplyr)
library(reshape)
library(ggplot2)
library(scales)
library(R.utils)

######### Helper objects for plotting #########

# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)


