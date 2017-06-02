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
sapply(unique(QC$PATIENT_ID), function(x){
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


######### VAFs of recurrent FF and FFPE variants #########

### Are there any recurrent variants? What are their VAFs?

# Collect all tier1-3 variants from all samples
allvars <- lapply(unique(QC$PATIENT_ID), function(x){
            # Load variants and tag them with PATIENT_ID
            filePath <- paste0("./Data/SNV/", x, "_SNVs_.tsv")
            var <- read.table(filePath, sep = "\t", header = T)
            var$PATIENT_ID <- as.character(x)
            return(var)
            })
allvars <- bind_rows(allvars)

# Make a table with variant KEY and how many times was observed in FFPE
allvars_ffpe <- allvars %>% filter(SAMPLE_TYPE == "FFPE")
recurr_ffpe <- lapply(unique(allvars_ffpe$KEY), function(x){
                 data.frame(KEY = x, NUM_OBS_FFPE = dim(allvars_ffpe[allvars_ffpe$KEY == x,])[1])
                  }
                  )
recurr_ffpe <- bind_rows(recurr_ffpe)
# Check for recurrent observations
table(recurr_ffpe$NUM_OBS_FFPE)
sum(recurr_ffpe$NUM_OBS_FFPE >1)  # 1566
sum(recurr_ffpe$NUM_OBS_FFPE >3)  # >5% samples, 415


# Make a table with variant KEY and how many times was observed in FF
allvars_ff <- allvars %>% filter(SAMPLE_TYPE == "FF")
recurr_ff <- lapply(unique(allvars_ff$KEY), function(x){
  data.frame(KEY = x, NUM_OBS_FF = dim(allvars_ff[allvars_ff$KEY == x,])[1])
}
)
recurr_ff <- bind_rows(recurr_ff)
# Check for recurrent observations
table(recurr_ff$NUM_OBS_FF)
sum(recurr_ff$NUM_OBS_FF >1)  # 878
sum(recurr_ff$NUM_OBS_FF >3)  # >5% samples, 70

# Check if FF recurrent variants are a subset of FFPE
recurr_ff_keys <- recurr_ff[recurr_ff$NUM_OBS_FF >3,]$KEY
recurr_ffpe_keys <- recurr_ffpe[recurr_ffpe$NUM_OBS_FFPE >3,]$KEY
sum(recurr_ff_keys %in% recurr_ffpe[recurr_ffpe$NUM_OBS_FFPE >3,]$KEY)  # 45 recurrent in both FF and FFPE

# Make a table of recurrent variants, flag those recurrent in FF, FFPE or both (recurrent in >5%)
recurr_subset <- allvars %>% filter(KEY %in% recurr_ff_keys | KEY %in% recurr_ffpe_keys)
recurr_subset$RECURR_TYPE <- ""
recurr_subset[(recurr_subset$KEY %in% recurr_ff_keys) & (recurr_subset$KEY %in% recurr_ffpe_keys), ]$RECURR_TYPE <- "BOTH"
recurr_subset[(recurr_subset$RECURR_TYPE != "BOTH") & (recurr_subset$KEY %in% recurr_ffpe_keys), ]$RECURR_TYPE <- "FFPE"
recurr_subset[(recurr_subset$RECURR_TYPE != "BOTH") & (recurr_subset$KEY %in% recurr_ff_keys), ]$RECURR_TYPE <- "FF"

# Write out the table with all observations of recurrent variants
write.table(recurr_subset, file = "./Data/SNV/var_recurrent_FF_FFPE_BOTH.tsv", sep = "\t", row.names = F, quote = F)

# Write out the table with basic info of unique recurrent variants
recurr <- recurr_subset %>% select(-(PATIENT_ID), -(SAMPLE_TYPE), -(VAF))
sum(duplicated(recurr))
recurr <- recurr[!duplicated(recurr),]
write.table(recurr, file = "./Data/SNV/var_recurrent_unique.tsv", sep = "\t", row.names = F, quote = F)


# Plot VAF distribution of recurrent variants (colored by whether are recurrent only in FF, FFPE or BOTH)
pdf(file = paste0("./Plots/Purity_VAF/VAF_recurr.pdf"))
print(ggplot(recurr_subset, aes(x = VAF, fill = factor(RECURR_TYPE))) + 
        geom_histogram(alpha = 0.5, position = "identity") + 
        scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) +
        theme(legend.position = "top") +
        blank)
dev.off()

# Exploratory analysis of recurrent variants
table(recurr$tier, recurr$class)
table(recurr$tier, recurr$RECURR_TYPE)

# Recurrent variants in Tier 1
recurr %>% filter(tier == "TIER1")
# Two recurrent KRAS mutations appearing mostly in colorectal samples
allvars %>% filter(KEY == "12_25245350_C_A")
allvars %>% filter(KEY == "12_25245347_C_T")
# One likely false positive in Tier1 (splice site) - problem with multi-allelic splits? (no VAFs listed)
allvars %>% filter(KEY == "2_47414419_GT_G")
# Tumor types where this variant is observed
QC %>% filter(PATIENT_ID %in% unique(allvars %>% filter(KEY == "2_47414419_GT_G") %>% .$PATIENT_ID)) %>% select(PATIENT_ID, TumorType)

# Recurrent variants in Tier 2
recurr %>% filter(tier == "TIER2")
# One that is recurrent only in FF
allvars %>% filter(KEY == "12_22625379_C_CG")
# Tumor types
QC %>% filter(PATIENT_ID %in% unique(allvars %>% filter(KEY == "12_22625379_C_CG") %>% .$PATIENT_ID)) %>% select(PATIENT_ID, TumorType)
# One recurrent in BOTH
allvars %>% filter(KEY == "4_53453080_CAG_C")





#########  Variant overlap plots ######### 

### Overlap of variants with no filtering

### Overlap of variants with recurrent variants removed 




#########  QC metrics correlation plots #########

### No variant filtering

### Recurrent variants removed 



