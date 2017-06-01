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


######### Load QC data #########

QC <- read.csv("./Data/QC_all62_FFPE_trios_clean_withPaths.csv")
QC <- QC %>% select(-(X))   # remove rownames written as 1st column



######### Tumor purity and allele frequency distributions #########


# For each patient, draw a distribution of tier1-3 allele frequencies (FF, FFPE groups) and indicate Illumina purity estimate
sapply(QC$PATIENT_ID, function(x){
  # Load variants
  filePath <- paste0("./Data/SNV/", x, "_SNVs_.tsv")
  var <- read.table(filePath, sep = "\t", header = T)
  # Draw histogram with lines at Illumina purity estimate, save plot
  pdf(file = paste0("./Plots/Purity_VAF/", x, "_VAF_purity.pdf"))
  print(ggplot(var, aes(x = VAF, fill = factor(SAMPLE_TYPE))) + 
      geom_histogram(alpha = 0.5, position = "identity") + 
      geom_vline(xintercept = ((QC %>% filter(PATIENT_ID == x, SAMPLE_TYPE == "FF") %>% .$TUMOUR_PURITY)/100), col = "coral", size = 1) + 
      geom_vline(xintercept = ((QC %>% filter(PATIENT_ID == x, SAMPLE_TYPE == "FFPE") %>% .$TUMOUR_PURITY)/100), col = "cyan3", size = 1) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) +
      theme(legend.position = "top") +
      blank)
  dev.off()
})






# Plot tumor purities




