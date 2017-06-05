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


# For each patient, draw a distribution of tier1-3 allele frequencies, VAF>0, (FF, FFPE groups) and indicate Illumina purity estimate
sapply(unique(QC$PATIENT_ID), function(x){
  # Load variants
  filePath <- paste0("./Data/SNV/", x, "_SNVs_0.tsv")
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
            filePath <- paste0("./Data/SNV/", x, "_SNVs_0.tsv")
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
sum(recurr_ffpe$NUM_OBS_FFPE >1)  # 1565
sum(recurr_ffpe$NUM_OBS_FFPE >3)  # >5% samples, 414


# Make a table with variant KEY and how many times was observed in FF
allvars_ff <- allvars %>% filter(SAMPLE_TYPE == "FF")
recurr_ff <- lapply(unique(allvars_ff$KEY), function(x){
  data.frame(KEY = x, NUM_OBS_FF = dim(allvars_ff[allvars_ff$KEY == x,])[1])
}
)
recurr_ff <- bind_rows(recurr_ff)
# Check for recurrent observations
table(recurr_ff$NUM_OBS_FF)
sum(recurr_ff$NUM_OBS_FF >1)  # 877
sum(recurr_ff$NUM_OBS_FF >3)  # >5% samples, 69

# Check if FF recurrent variants are a subset of FFPE
recurr_ff_keys <- recurr_ff[recurr_ff$NUM_OBS_FF >3,]$KEY
recurr_ffpe_keys <- recurr_ffpe[recurr_ffpe$NUM_OBS_FFPE >3,]$KEY
sum(recurr_ff_keys %in% recurr_ffpe[recurr_ffpe$NUM_OBS_FFPE >3,]$KEY)  # 44 recurrent in both FF and FFPE

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
# Add the mean VAF for each recurrent variant
recurr$MEAN_VAF <- sapply(recurr$KEY, function(x){
  mean(recurr_subset[recurr_subset$KEY == x,]$VAF)
})
# Plot mean VAF of recurrent variants
pdf(file = paste0("./Plots/Purity_VAF/VAF_recurr_byTier.pdf"))
print(ggplot(recurr, aes(x=RECURR_TYPE, y=MEAN_VAF, fill = factor(tier))) + 
        geom_boxplot() + 
        labs(x = "Sample type", y = "Mean VAF") + 
        blank)
dev.off()
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
# Two recurrent KRAS mutations appearing mostly in colorectal samples, high VAFs, in both FF and FFPE
allvars %>% filter(KEY == "12_25245350_C_A")
allvars %>% filter(KEY == "12_25245347_C_T")
# Third recurrent
allvars %>% filter(KEY == "3_179234297_A_G")  # known PIK3CA mutation in breast cancer, high VAFs

# Tumor types where this variant is observed
QC %>% filter(PATIENT_ID %in% unique(allvars %>% filter(KEY == "12_25245350_C_A") %>% .$PATIENT_ID)) %>% select(PATIENT_ID, TumorType)
QC %>% filter(PATIENT_ID %in% unique(allvars %>% filter(KEY == "12_25245347_C_T") %>% .$PATIENT_ID)) %>% select(PATIENT_ID, TumorType)
QC %>% filter(PATIENT_ID %in% unique(allvars %>% filter(KEY == "3_179234297_A_G") %>% .$PATIENT_ID)) %>% select(PATIENT_ID, TumorType)


# Recurrent variants in Tier 2
recurr %>% filter(tier == "TIER2")
# One that is recurrent only in FF
allvars %>% filter(KEY == "12_22625379_C_CG")  # both in FF and FFPE, low VAF, repeats region
# Tumor types
QC %>% filter(PATIENT_ID %in% unique(allvars %>% filter(KEY == "12_22625379_C_CG") %>% .$PATIENT_ID)) %>% select(PATIENT_ID, TumorType)
# One recurrent in BOTH
allvars %>% filter(KEY == "4_53453080_CAG_C")  # mostly in FFPE, within simple repeats, likely false positive


# Plot all variants (recurrent and non-recurrent) mean VAFs in FF and FFPE 
pdf(file = paste0("./Plots/Purity_VAF/VAF_all_byTier.pdf"))
print(ggplot(allvars, aes(x=SAMPLE_TYPE, y=VAF, fill = factor(tier))) + 
        geom_boxplot() + 
        labs(x = "Sample type", y = "VAF") + 
        blank)
dev.off()





#########  Variant overlap plots ######### 

### Overlap of variants with no filtering

# Load SNV summary table for VAF > 0 
SNV_summary <- read.csv("./Data/SNV/SNV_summary_62trios_0_2017-06-05.csv")

# Add QC details to the SNV summary table
SNV_summary$PATIENT_ID <- as.character(SNV_summary$PATIENT_ID)
QC_table <- QC %>% filter(SAMPLE_TYPE == "FFPE") %>% dplyr::select(PATIENT_ID, CENTER_CODE, TumorType, SAMPLE_WELL_ID, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER, DEAMINATION_MISMATCHES_PER)
names(QC_table)[4:12] <- paste0("FFPE_", names(QC_table)[4:12])
QC_table <- left_join(QC_table, (QC %>% filter(SAMPLE_TYPE == "FF") %>% dplyr::select(PATIENT_ID, SAMPLE_WELL_ID, LIBRARY_TYPE, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER, DEAMINATION_MISMATCHES_PER)), by = "PATIENT_ID")
names(QC_table)[13:22] <- paste0("FF_", names(QC_table)[13:22])
QC_table$PATIENT_ID <- as.character(QC_table$PATIENT_ID)
SNV_summary <- left_join(SNV_summary, QC_table, by = "PATIENT_ID")

### Summary plots
# Distribution of RECALL and PRECISION
hist(SNV_summary$RECALL)
hist(SNV_summary$PRECISION)

### Overlap plot
# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
SNV_summary_m <- as.data.frame(t(SNV_summary %>% arrange(PRECISION) %>% dplyr::select(PATIENT_ID, OVERLAP, FF_UNIQ, FFPE_UNIQ)))
names(SNV_summary_m) <- as.matrix(SNV_summary_m[1,])
SNV_summary_m <- SNV_summary_m[2:4,]
SNV_summary_m <- melt(cbind(SNV_summary_m, ind = rownames(SNV_summary_m)), id.vars = c('ind'))
SNV_summary_m$value <- as.numeric(SNV_summary_m$value)
# Overap plot (needs package "scales")
pdf(file = paste0("./Plots/Concordance/Overlap_all.pdf"))
print(ggplot(SNV_summary_m, aes(x = variable, y = value, fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format()) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + 
  theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + 
  blank)
dev.off()


### Overlap of variants with FFPE recurrent variants, mean VAF < 0.1 removed 

# Load SNV summary table for VAF > 0, recurrent variants with VAF < 0.1 removed
SNV_summary_R <- read.csv("./Data/SNV/SNV_summary_62trios_noRecurr_0_2017-06-05.csv")

# Add QC details to the SNV summary table
SNV_summary_R$PATIENT_ID <- as.character(SNV_summary_R$PATIENT_ID)
SNV_summary_R <- left_join(SNV_summary_R, QC_table, by = "PATIENT_ID")

### Summary plots
# Distribution of RECALL and PRECISION
hist(SNV_summary_R$RECALL)
hist(SNV_summary_R$PRECISION)

### Overlap plot
# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
SNV_summary_m <- as.data.frame(t(SNV_summary_R %>% arrange(PRECISION) %>% dplyr::select(PATIENT_ID, OVERLAP, FF_UNIQ, FFPE_UNIQ)))
names(SNV_summary_m) <- as.matrix(SNV_summary_m[1,])
SNV_summary_m <- SNV_summary_m[2:4,]
SNV_summary_m <- melt(cbind(SNV_summary_m, ind = rownames(SNV_summary_m)), id.vars = c('ind'))
SNV_summary_m$value <- as.numeric(SNV_summary_m$value)
# Overap plot (needs package "scales")
pdf(file = paste0("./Plots/Concordance/Overlap_noRecurr.pdf"))
print(ggplot(SNV_summary_m, aes(x = variable, y = value, fill = ind)) + 
        geom_bar(position = "fill",stat = "identity") + 
        scale_y_continuous(labels = percent_format()) + 
        theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + 
        theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + 
        blank)
dev.off()


### Overlap with all VAF > 0.1

# Load SNV summary table for VAF > 0, recurrent variants with VAF < 0.1 removed
SNV_summary_R_10 <- read.csv("./Data/SNV/SNV_summary_62trios_noRecurr_0.1_2017-06-05.csv")

# Add QC details to the SNV summary table
SNV_summary_R_10$PATIENT_ID <- as.character(SNV_summary_R_10$PATIENT_ID)
SNV_summary_R_10 <- left_join(SNV_summary_R_10, QC_table, by = "PATIENT_ID")

### Summary plots
# Distribution of RECALL and PRECISION
hist(SNV_summary_R_10$RECALL)
hist(SNV_summary_R_10$PRECISION)

### Overlap plot
# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
SNV_summary_m <- as.data.frame(t(SNV_summary_R_10 %>% arrange(PRECISION) %>% dplyr::select(PATIENT_ID, OVERLAP, FF_UNIQ, FFPE_UNIQ)))
names(SNV_summary_m) <- as.matrix(SNV_summary_m[1,])
SNV_summary_m <- SNV_summary_m[2:4,]
SNV_summary_m <- melt(cbind(SNV_summary_m, ind = rownames(SNV_summary_m)), id.vars = c('ind'))
SNV_summary_m$value <- as.numeric(SNV_summary_m$value)
# Overap plot (needs package "scales")
pdf(file = paste0("./Plots/Concordance/Overlap_noRecurr_10.pdf"))
print(ggplot(SNV_summary_m, aes(x = variable, y = value, fill = ind)) + 
        geom_bar(position = "fill",stat = "identity") + 
        scale_y_continuous(labels = percent_format()) + 
        theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + 
        theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + 
        blank)
dev.off()




### Overlap with all VAF > 0.2

# Load SNV summary table for VAF > 0, recurrent variants with VAF < 0.1 removed
SNV_summary_R_20 <- read.csv("./Data/SNV/SNV_summary_62trios_noRecurr_0.2_2017-06-05.csv")

# Add QC details to the SNV summary table
SNV_summary_R_20$PATIENT_ID <- as.character(SNV_summary_R_20$PATIENT_ID)
SNV_summary_R_20 <- left_join(SNV_summary_R_20, QC_table, by = "PATIENT_ID")

### Summary plots
# Distribution of RECALL and PRECISION
hist(SNV_summary_R_20$RECALL)
hist(SNV_summary_R_20$PRECISION)

### Overlap plot
# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
SNV_summary_m <- as.data.frame(t(SNV_summary_R_10 %>% arrange(PRECISION) %>% dplyr::select(PATIENT_ID, OVERLAP, FF_UNIQ, FFPE_UNIQ)))
names(SNV_summary_m) <- as.matrix(SNV_summary_m[1,])
SNV_summary_m <- SNV_summary_m[2:4,]
SNV_summary_m <- melt(cbind(SNV_summary_m, ind = rownames(SNV_summary_m)), id.vars = c('ind'))
SNV_summary_m$value <- as.numeric(SNV_summary_m$value)
# Overap plot (needs package "scales")
pdf(file = paste0("./Plots/Concordance/Overlap_noRecurr_20.pdf"))
print(ggplot(SNV_summary_m, aes(x = variable, y = value, fill = ind)) + 
        geom_bar(position = "fill",stat = "identity") + 
        scale_y_continuous(labels = percent_format()) + 
        theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + 
        theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + 
        blank)
dev.off()



################ Boxplots of recall vs precission ################

# Gather data into one table
SNV_summary$GROUP <- "ALL"
SNV_summary_R$GROUP <- "NO_RECURRENT"
SNV_summary_R_10$GROUP <- "NO_RECURRENT_VAF_10"
SNV_summary_R_20$GROUP <- "NO_RECURRENT_VAF_20"

SNV_summary_all <- rbind(SNV_summary, SNV_summary_R, SNV_summary_R_10, SNV_summary_R_20)

# Plot recall
pdf(file = paste0("./Plots/Concordance/Recall_box.pdf"))
print(ggplot(SNV_summary_all, aes(x=GROUP, y=RECALL, fill = factor(GROUP))) + 
        geom_boxplot() + 
        labs(x="", y="RECALL") +
        theme(axis.text.x=element_blank()) + 
        blank)
dev.off()

# Plot precission
pdf(file = paste0("./Plots/Concordance/Precission_box.pdf"))
print(ggplot(SNV_summary_all, aes(x=GROUP, y=PRECISION, fill = factor(GROUP))) + 
        geom_boxplot() + 
        labs(x="", y="PRECISION") +
        theme(axis.text.x=element_blank()) + 
        blank)
dev.off()


################ FF vs FFPE VAF plot ################ 

# Load all somatic variants in FF and FFPE for BRC patient 200000306 (FF = LP2000407-DNA_A01, FFPE = LP2000418-DNA_A01)
ff <- read.table("./Data/SNV/LP2000407-DNA_A01-info.txt", sep = "\t")
ffpe <- read.table("./Data/SNV/LP2000418-DNA_A01-info.txt", sep = "\t")
names(ff) <- c("chr", "pos", "ref", "alt", "VAF")
names(ffpe) <- c("chr", "pos", "ref", "alt", "VAF")
ff$GROUP <- "FF"
ffpe$GROUP <- "FFPE"

# Create keys to compare variants
ff$chr <- sub("chr", "", ff$chr)
ffpe$chr <- sub("chr", "", ffpe$chr)
ff$KEY <- sapply(1:dim(ff)[1], function(x){
  paste(ff$chr[x], ff$pos[x], ff$ref[x], ff$alt[x], sep = "_")
})
ffpe$KEY <- sapply(1:dim(ffpe)[1], function(x){
  paste(ffpe$chr[x], ffpe$pos[x], ffpe$ref[x], ffpe$alt[x], sep = "_")
})

# Merge FF and FFPE
ff_ffpe <- rbind(ff, ffpe)

# Plot VAF distribution for the whole genome variants
# Draw histogram with lines at Illumina purity estimate, save plot
pdf(file = paste0("./Plots/Purity_VAF/200000306_VAF_purity_allVariants.pdf"))
print(ggplot(ff_ffpe, aes(x = VAF, fill = factor(GROUP))) + 
        geom_histogram(alpha = 0.5, position = "identity") + 
        geom_vline(xintercept = ((QC %>% filter(PATIENT_ID == "200000306", SAMPLE_TYPE == "FF") %>% .$TUMOUR_PURITY)/100), col = "coral", size = 1) + 
        geom_vline(xintercept = ((QC %>% filter(PATIENT_ID == "200000306", SAMPLE_TYPE == "FFPE") %>% .$TUMOUR_PURITY)/100), col = "cyan3", size = 1) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) +
        theme(legend.position = "top") +
        blank)
dev.off()


# Summarize by KEY (set NAs to zero); KEY in rows, VAF_FF and VAF_FFPE in columns
names(ff)[5] <- "VAF_FF"
names(ffpe)[5] <- "VAF_FFPE"
vaf_table <- full_join((ff %>% select(KEY, VAF_FF)), (ffpe %>% select(KEY, VAF_FFPE)), by = "KEY")
vaf_table[is.na(vaf_table$VAF_FF),]$VAF_FF <- 0
vaf_table[is.na(vaf_table$VAF_FFPE),]$VAF_FFPE <- 0


# Plot VAF, FF vs FFPE
pdf(file = paste0("./Plots/Purity_VAF/200000306_VAF_allVariants.pdf"))
print(ggplot(vaf_table, aes(x=VAF_FF, y=VAF_FFPE)) +
  geom_point(size = 3, alpha = 0.02, col = "violet") +
  blank)
dev.off()


### Overlap of variants BY TIER

tier1 <- unique(allvars %>% filter(PATIENT_ID == "200000306", tier == "TIER1") %>% .$KEY)
tier2 <- unique(allvars %>% filter(PATIENT_ID == "200000306", tier == "TIER2") %>% .$KEY)
tier3 <- unique(allvars %>% filter(PATIENT_ID == "200000306", tier == "TIER3") %>% .$KEY)

vaf_table$TIER <- "NONE"
vaf_table[vaf_table$KEY %in% tier1,]$TIER <- "TIER1"
vaf_table[vaf_table$KEY %in% tier2,]$TIER <- "TIER2"
vaf_table[vaf_table$KEY %in% tier3,]$TIER <- "TIER3"

pdf(file = paste0("./Plots/Purity_VAF/200000306_VAF_tier123.pdf"))
print(ggplot(vaf_table %>% filter(TIER != "NONE"), aes(x=VAF_FF, y=VAF_FFPE, col = factor(TIER))) +
        geom_point(size = 3, alpha = 0.5) +
        blank)
dev.off()



#########  QC metrics correlation plots #########

# Mark BRC pilot samples
table(SNV_summary$CENTER_CODE, exclude = NULL)
SNV_summary$CENTER_CODE <- as.character(SNV_summary$CENTER_CODE)
SNV_summary[is.na(SNV_summary$CENTER_CODE),]$CENTER_CODE <- "BRC"

SNV_summary_R_10$CENTER_CODE <- as.character(SNV_summary_R_10$CENTER_CODE)
SNV_summary_R_10[is.na(SNV_summary_R_10$CENTER_CODE),]$CENTER_CODE <- "BRC"


### No variant filtering

# RECALL of FF by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, etc
ggplot(SNV_summary, aes(x = RECALL, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_AT_DROP, method = "spearman")  # 0.01094204

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # 0.06393162

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_CHIMERIC_PER, method = "spearman")  # -0.08038311

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.07213675

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # -0.1070085

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(SNV_summary$RECALL, SNV_summary$FFPE_TUMOUR_PURITY, method = "spearman")  # -0.3930165


# PRECISION of FF by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, etc
ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_AT_DROP, method = "spearman")  # -0.3135579

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # -0.3996581

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_CHIMERIC_PER, method = "spearman")  # -0.04241492

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.2205128

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # 0.04820513

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Precision", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(SNV_summary$PRECISION, SNV_summary$FFPE_TUMOUR_PURITY, method = "spearman")  # -0.4334136




### FFPE recurrent variants, mean VAF < 0.1 removed

# RECALL of FF by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, etc
ggplot(SNV_summary_R_10, aes(x = RECALL, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$RECALL, SNV_summary_R_10$FFPE_AT_DROP, method = "spearman")  # 0.01094204

ggplot(SNV_summary_R_10, aes(x = RECALL, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$RECALL, SNV_summary_R_10$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # 0.06393162

ggplot(SNV_summary_R_10, aes(x = RECALL, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$RECALL, SNV_summary_R_10$FFPE_CHIMERIC_PER, method = "spearman")  # -0.08038311

ggplot(SNV_summary_R_10, aes(x = RECALL, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$RECALL, SNV_summary_R_10$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.07213675

ggplot(SNV_summary_R_10, aes(x = RECALL, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$RECALL, SNV_summary_R_10$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # -0.1070085

ggplot(SNV_summary_R_10, aes(x = RECALL, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(SNV_summary_R_10$RECALL, SNV_summary_R_10$FFPE_TUMOUR_PURITY, method = "spearman")  # -0.3930165




# PRECISION of FF by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, etc
ggplot(SNV_summary_R_10, aes(x = PRECISION, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$PRECISION, SNV_summary_R_10$FFPE_AT_DROP, method = "spearman")  # -0.3135579

ggplot(SNV_summary_R_10, aes(x = PRECISION, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$PRECISION, SNV_summary_R_10$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # -0.3996581

ggplot(SNV_summary_R_10, aes(x = PRECISION, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$PRECISION, SNV_summary_R_10$FFPE_CHIMERIC_PER, method = "spearman")  # -0.04241492

ggplot(SNV_summary_R_10, aes(x = PRECISION, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$PRECISION, SNV_summary_R_10$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.2205128

ggplot(SNV_summary_R_10, aes(x = PRECISION, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary_R_10$PRECISION, SNV_summary_R_10$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # 0.04820513

ggplot(SNV_summary_R_10, aes(x = PRECISION, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Precision", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(SNV_summary_R_10$PRECISION, SNV_summary_R_10$FFPE_TUMOUR_PURITY, method = "spearman")  # -0.4334136




############### Analysing well performing samples ############### 

summary(SNV_summary_R_10 %>% filter(RECALL > 0.75 & PRECISION > 0.75))
SNV_summary_R_10 %>% filter(RECALL > 0.75 & PRECISION > 0.75) %>% select(PATIENT_ID, TumorType, FF_TUMOUR_PURITY, FFPE_TUMOUR_PURITY, FFPE_AT_DROP, FFPE_COVERAGE_HOMOGENEITY)

# Recall and precision colored by various FFPE metrics
ggplot(SNV_summary_R_10, aes(x = RECALL, y = PRECISION, col = factor(RECALL > 0.75 & PRECISION > 0.75))) + geom_jitter() + labs(x = "Recall", y = "Precision")  + theme(legend.title=element_blank()) + regr_line
ggplot(SNV_summary_R_10, aes(x = RECALL, y = PRECISION, col = TumorType)) + geom_jitter() + labs(x = "Recall", y = "Precision")  + theme(legend.title=element_blank()) + regr_line
ggplot(SNV_summary_R_10, aes(x = RECALL, y = PRECISION, col = FFPE_AT_DROP)) + geom_jitter() + labs(x = "Recall", y = "Precision")  + theme(legend.title=element_blank()) + regr_line
ggplot(SNV_summary_R_10, aes(x = RECALL, y = PRECISION, col = FFPE_COVERAGE_HOMOGENEITY)) + geom_jitter() + labs(x = "Recall", y = "Precision")  + theme(legend.title=element_blank()) + regr_line
ggplot(SNV_summary_R_10, aes(x = RECALL, y = PRECISION, col = FFPE_TUMOUR_PURITY)) + geom_jitter() + labs(x = "Recall", y = "Precision")  + theme(legend.title=element_blank()) + regr_line
ggplot(SNV_summary_R_10, aes(x = RECALL, y = PRECISION, col = FF_TUMOUR_PURITY)) + geom_jitter() + labs(x = "Recall", y = "Precision")  + theme(legend.title=element_blank()) + regr_line
ggplot(SNV_summary_R_10, aes(x = RECALL, y = PRECISION, col = FFPE_DEAMINATION_MISMATCHES_PER)) + geom_jitter() + labs(x = "Recall", y = "Precision")  + theme(legend.title=element_blank()) + regr_line


# Samples performing particularly badly
SNV_summary_R_10 %>% filter(RECALL < 0.3 & PRECISION < 0.3) %>% select(PATIENT_ID, TumorType, FF_TUMOUR_PURITY, FFPE_TUMOUR_PURITY, FFPE_AT_DROP, FFPE_COVERAGE_HOMOGENEITY)




