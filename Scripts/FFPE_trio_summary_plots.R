# FFPE trio analysis
# Martina Mijuskovic
# May 2017

######### Helper objects for plotting #########

# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)



######### QC summary plots ######### 


######### Calculate SNV overlaps ######### 


######### Calculate SV overlaps ######### 


######### Calculate CNV overlaps ######### 


######### Overlap plots ######### 


######### Correlation plots ######### 

