# Martina Mijuskovic
# FFPE project
# Compares SNV calls between FFPE and FF samples from FFPE trios (HPC)

# NOTE: Installing R packages on HPC: use lib = "~/R/x86_64-pc-linux-gnu-library/3.3"

library(dplyr)
library(jsonlite)
library(data.table)

# Working directory on the HPC
setwd("/home/mmijuskovic/FFPE_trio_analysis/SNV_trio_comparison")

# Today's date
today <- Sys.Date()



############# Get SNV VCF and JSON paths ############# 

# Load the manifest (HPC)
QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE_trio_analysis/QC_all62_FFPE_trios_clean.csv")  # all 62 trios
# QC_portal_trios <- read.csv("./Data/QC_portal_trios_final.csv") # local
# QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE_trio_analysis/QC_all62_FFPE_trios_clean_withPaths.csv")  # already subset, containing VCF & JSON paths
# QC_portal_trios$json_path <- as.character(QC_portal_trios$json_path)  # It reads in as factor, causing errors downstream

# Subset for FF and FFPE samples
QC_portal_trios <- QC_portal_trios %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE"))

# Get VCF paths (HPC)   ------ original Illumina SNVs VCF! --- some folders contain two, one ending in "PASS.duprem.atomic.left.split.somatic.vcf.gz"
paths <- unlist(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.somatic.vcf.gz", sep = " ")
  system(command, intern = T)
}))
paths <- paths[!grepl("PASS.duprem.atomic.left.split.somatic.vcf.gz", paths)]
QC_portal_trios$SNV_VCF_path <- as.character(paths)

# Get SNV tiering json paths (HPC)
system("touch /home/mmijuskovic/FFPE/FFPE_trios_timestamp -d 2017-03-29") # timestamp to help select only new tiering json files (produced after March 29th 2017)
QC_portal_trios$json_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", sub("/genomes", "/genomes/analysis", x), "-newer /home/mmijuskovic/FFPE/FFPE_trios_timestamp -iname *tiering.json", sep = " ")
  system(command, intern = T)
}))

# Write the QC portal table with paths
write.csv(QC_portal_trios, file = "/home/mmijuskovic/FFPE_trio_analysis/QC_all62_FFPE_trios_clean_withPaths.csv")


############# Calculate SNV concordance ############# 

### Function that extracts variant information from json files and calculates concordance between FF and FFPE
compareSNV <- function(patientID, var_freq = NULL){
  
  # Get patientID
  ID <- patientID
  
  # Get json paths
  ff_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FF") %>% .$json_path
  ffpe_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FFPE") %>% .$json_path
  
  ### Read json files: get chromosome, position, VAF, tier, etc
  
  # FF
  ff_json <- fromJSON(ff_path, flatten = T)
  ff <- ff_json %>% dplyr::select(
    reportedVariantCancer.chromosome, reportedVariantCancer.position, reportedVariantCancer.reference, reportedVariantCancer.alternate, 
    reportedVariantCancer.VAF, somaticOrGermline)
  # Reduce to somatic, remove that column
  ff <- ff %>% filter(somaticOrGermline == "somatic") %>% dplyr::select(-(somaticOrGermline))
  ff$tier <- sapply(1:dim(ff)[1], function(x){
    ff_json$reportedVariantCancer.reportEvents[[x]]$tier})
  ff$class <- sapply(1:dim(ff)[1], function(x){
    ff_json$reportedVariantCancer.reportEvents[[x]]$soNames})
  names(ff) <- c("chr", "pos", "ref", "alt", "VAF", "tier", "class")
  ff$KEY <- sapply(1:dim(ff)[1], function(x){
    paste(ff$chr[x], ff$pos[x], ff$ref[x], ff$alt[x], sep = "_")
  })

  # FFPE
  ffpe_json <- fromJSON(ffpe_path, flatten = T)
  ffpe <- ffpe_json %>% dplyr::select(
    reportedVariantCancer.chromosome, reportedVariantCancer.position, reportedVariantCancer.reference, reportedVariantCancer.alternate, 
    reportedVariantCancer.VAF, somaticOrGermline)
  # Reduce to somatic, remove that column
  ffpe <- ffpe %>% filter(somaticOrGermline == "somatic") %>% dplyr::select(-(somaticOrGermline))
  ffpe$tier <- sapply(1:dim(ffpe)[1], function(x){
    ffpe_json$reportedVariantCancer.reportEvents[[x]]$tier})
  ffpe$class <- sapply(1:dim(ffpe)[1], function(x){
    ffpe_json$reportedVariantCancer.reportEvents[[x]]$soNames})
  names(ffpe) <- c("chr", "pos", "ref", "alt", "VAF", "tier", "class")
  ffpe$KEY <- sapply(1:dim(ffpe)[1], function(x){
    paste(ffpe$chr[x], ffpe$pos[x], ffpe$ref[x], ffpe$alt[x], sep = "_")
  }) 
 
  # # Sanity check (!!! found duplicates - traced the error to the tiering pipeline)
  # sum(duplicated(ff))  #33  # after OpenCGA fix -> 1
  # sum(duplicated(ff$KEY))  #37  # after OpenCGA fix -> 1
  # # Checking the leftover duplicate
  # dup_key <- ff[duplicated(ff$KEY),]$KEY
  # ff %>% filter(KEY == dup_key)

  # Deduplicate
  ff <- ff[!duplicated(ff),]
  ffpe <- ffpe[!duplicated(ffpe),]
  
  # # Check duplicate keys
  # sum(duplicated(ff$KEY))
  # sum(duplicated(ffpe$KEY))
  # ff %>% filter(KEY %in% (ff[duplicated(ff$KEY),]$KEY))
  
  # Deduplicate variants with same keys
  ff <- ff[(!duplicated(ff$KEY)),]
  ffpe <- ffpe[(!duplicated(ffpe$KEY)),]
  
  # Subset the variants by frequency if frequency argument is provided
  if (!is.null(var_freq)) {
    ff <- ff %>% filter(VAF > var_freq)
    ffpe <- ffpe %>% filter(VAF > var_freq)
  }
  
  # Write the table with filtered variants
  if (dim(ff)[1] != 0) {
    ff$SAMPLE_TYPE <- "FF"}
  if (dim(ffpe)[1] != 0) {
    ffpe$SAMPLE_TYPE <- "FFPE"
  }
  all <- rbind(ff, ffpe)
  all$class <- as.character(all$class)
  write.table(all, file = paste0(ID, "_SNVs_", var_freq, ".tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
  
  # Summary table with concordance
  result <- data.frame(FF_TOTAL = dim(ff)[1],
                       FFPE_TOTAL = dim(ffpe)[1],
                       OVERLAP = sum(ff$KEY %in% ffpe$KEY), # bug
                       FF_UNIQ = ((dim(ff)[1]) - sum(ff$KEY %in% ffpe$KEY)),
                       FFPE_UNIQ = ((dim(ffpe)[1]) - sum(ff$KEY %in% ffpe$KEY)),
                       RECALL = ((sum(ff$KEY %in% ffpe$KEY))/(dim(ff)[1])),
                       PRECISION = (sum(ff$KEY %in% ffpe$KEY))/(dim(ffpe)[1]))

 return(result)
}

# Get all patient IDs
patientIDs <- as.character(unique(QC_portal_trios$PATIENT_ID))

# Run SNV comparison for each patient ID (all allele frequencies)
SNV_summary <- lapply(patientIDs, compareSNV)
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs

# Write out the resulting table
write.csv(SNV_summary, file = paste0("SNV_summary_62trios_allFreq_", today, ".csv"), row.names = F, quote = F)


# Run with different variant frequency thresholds
freq <- seq(from = 0, to = 0.3, by = 0.05)
lapply(freq, function(x){
  SNV_summary <- lapply(patientIDs, compareSNV, var_freq = x)
  SNV_summary <- bind_rows(SNV_summary)
  SNV_summary$PATIENT_ID <- patientIDs
  write.csv(SNV_summary, file = paste0("SNV_summary_62trios_", x, "_", today, ".csv"), row.names = F, quote = F)
})


# NOTE that file "allFreq" and "0" are slightly different
allFreq <- read.csv("./Data/SNV/SNV_summary_62trios_allFreq_2017-06-05.csv")
zeroFreq <- read.csv("./Data/SNV/SNV_summary_62trios_0_2017-06-05.csv")
allFreq == zeroFreq
# > allFreq %>% filter(PATIENT_ID == "200000312")
# FF_TOTAL FFPE_TOTAL OVERLAP FF_UNIQ FFPE_UNIQ    RECALL PRECISION PATIENT_ID
# 1       36        107      19      17        88 0.5277778 0.1775701  200000312
# > zeroFreq  %>% filter(PATIENT_ID == "200000312")
# FF_TOTAL FFPE_TOTAL OVERLAP FF_UNIQ FFPE_UNIQ    RECALL PRECISION PATIENT_ID
# 1       35        107      18      17        89 0.5142857 0.1682243  200000312
# > allFreq %>% filter(PATIENT_ID == "200000335")
# FF_TOTAL FFPE_TOTAL OVERLAP FF_UNIQ FFPE_UNIQ    RECALL PRECISION PATIENT_ID
# 1      102        276      66      36       210 0.6470588 0.2391304  200000335
# > zeroFreq  %>% filter(PATIENT_ID == "200000335")
# FF_TOTAL FFPE_TOTAL OVERLAP FF_UNIQ FFPE_UNIQ    RECALL PRECISION PATIENT_ID
# 1      102        275      66      36       209 0.6470588      0.24  200000335

# Finding the source of discrepancy in patient 200000312
all_200000312 <- fread("./Data/SNV/200000312_SNVs_.tsv", sep = "\t")
zero_200000312 <- read.table("./Data/SNV/200000312_SNVs_0.tsv", skip = 1, sep = "\t", colClasses = c("character", "integer", "character", "character", "numeric", rep("character",4)))
names(zero_200000312) <- names(all_200000312)
all_equal(all_200000312, zero_200000312 )  # "Different number of rows"
table(zero_200000312$SAMPLE_TYPE)
table(all_200000312$SAMPLE_TYPE)
# Finding the rows missing in FF (but present in FFPE)
all_200000312[all_200000312$SAMPLE_TYPE == "FF",]$KEY[!all_200000312[all_200000312$SAMPLE_TYPE == "FF",]$KEY %in% zero_200000312[zero_200000312$SAMPLE_TYPE == "FF",]$KEY]
all_200000312 %>% filter(KEY == "13_20082011_TA_T")
zero_200000312 %>% filter(KEY == "13_20082011_TA_T")

# Finding the source of discrepancy in patient 200000335
all_200000335 <- fread("./Data/SNV/200000335_SNVs_.tsv", sep = "\t")
zero_200000335 <- read.table("./Data/SNV/200000335_SNVs_0.tsv", skip = 1, sep = "\t", colClasses = c("character", "integer", "character", "character", "numeric", rep("character",4)))
names(zero_200000335) <- names(all_200000335)
all_equal(all_200000335, zero_200000335)  # "Different number of rows"
table(zero_200000335$SAMPLE_TYPE)
table(all_200000335$SAMPLE_TYPE)
# Finding the rows missing in FFPE (but present in FF)
all_200000335[all_200000335$SAMPLE_TYPE == "FFPE",]$KEY[!all_200000335[all_200000335$SAMPLE_TYPE == "FFPE",]$KEY %in% zero_200000335[zero_200000335$SAMPLE_TYPE == "FFPE",]$KEY]
all_200000335 %>% filter(KEY == "1_21480128_T_C")  # extra variant in FFPE, VAF = 0
zero_200000335 %>% filter(KEY == "1_21480128_T_C")


# CONCLUSION: 
# Looking at the original and normalized VCF, it became clear that this is a Illumina's subtraction bug. 
# The indel variant has read depth reported in the "TAR" format field instead of "TIR" (another bug). 
# Since this is an indel, and TIR is zero, our calculated VAF is zero.
# The other variant has filtered depth of zero in germline and only 5 in somatic. It should't PASS.
# Issue needs to be reported to Illumina.
# I will only use files with VAF > 0 for FF/FFPE comparisons.



############# Calculate SNV concordance removing low VAF recurrent variants ############# 

# Load the table with recurrent variants and mean VAFs (created in SNV_plots.R)
#recurr <- read.table("./Data/SNV/var_recurrent_unique.tsv", sep = "\t", header = T) # local
recurr <- read.table("/home/mmijuskovic/FFPE_trio_analysis/SNV_trio_comparison/var_recurrent_unique.tsv", sep = "\t", header = T)

# Get KEYs of variants with VAF < 0.1
recurr_keys <- as.character(recurr %>% filter(MEAN_VAF < 0.1) %>% .$KEY)  # 240 variants

# Calculate concordance again, removing recurrent variants, no variant tables written
compareSNV_2 <- function(patientID, var_freq = NULL){
  
  # Get patientID
  ID <- patientID
  
  # Get json paths
  ff_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FF") %>% .$json_path
  ffpe_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FFPE") %>% .$json_path
  
  ### Read json files: get chromosome, position, VAF, tier, etc
  
  # FF
  ff_json <- fromJSON(ff_path, flatten = T)
  ff <- ff_json %>% dplyr::select(
    reportedVariantCancer.chromosome, reportedVariantCancer.position, reportedVariantCancer.reference, reportedVariantCancer.alternate, 
    reportedVariantCancer.VAF, somaticOrGermline)
  # Reduce to somatic, remove that column
  ff <- ff %>% filter(somaticOrGermline == "somatic") %>% dplyr::select(-(somaticOrGermline))
  ff$tier <- sapply(1:dim(ff)[1], function(x){
    ff_json$reportedVariantCancer.reportEvents[[x]]$tier})
  ff$class <- sapply(1:dim(ff)[1], function(x){
    ff_json$reportedVariantCancer.reportEvents[[x]]$soNames})
  names(ff) <- c("chr", "pos", "ref", "alt", "VAF", "tier", "class")
  ff$KEY <- sapply(1:dim(ff)[1], function(x){
    paste(ff$chr[x], ff$pos[x], ff$ref[x], ff$alt[x], sep = "_")
  })
  
  # FFPE
  ffpe_json <- fromJSON(ffpe_path, flatten = T)
  ffpe <- ffpe_json %>% dplyr::select(
    reportedVariantCancer.chromosome, reportedVariantCancer.position, reportedVariantCancer.reference, reportedVariantCancer.alternate, 
    reportedVariantCancer.VAF, somaticOrGermline)
  # Reduce to somatic, remove that column
  ffpe <- ffpe %>% filter(somaticOrGermline == "somatic") %>% dplyr::select(-(somaticOrGermline))
  ffpe$tier <- sapply(1:dim(ffpe)[1], function(x){
    ffpe_json$reportedVariantCancer.reportEvents[[x]]$tier})
  ffpe$class <- sapply(1:dim(ffpe)[1], function(x){
    ffpe_json$reportedVariantCancer.reportEvents[[x]]$soNames})
  names(ffpe) <- c("chr", "pos", "ref", "alt", "VAF", "tier", "class")
  ffpe$KEY <- sapply(1:dim(ffpe)[1], function(x){
    paste(ffpe$chr[x], ffpe$pos[x], ffpe$ref[x], ffpe$alt[x], sep = "_")
  }) 
  
  # # Sanity check (!!! found duplicates - traced the error to the tiering pipeline)
  # sum(duplicated(ff))  #33  # after OpenCGA fix -> 1
  # sum(duplicated(ff$KEY))  #37  # after OpenCGA fix -> 1
  # # Checking the leftover duplicate
  # dup_key <- ff[duplicated(ff$KEY),]$KEY
  # ff %>% filter(KEY == dup_key)
  
  # Deduplicate
  ff <- ff[!duplicated(ff),]
  ffpe <- ffpe[!duplicated(ffpe),]
  
  # # Check duplicate keys
  # sum(duplicated(ff$KEY))
  # sum(duplicated(ffpe$KEY))
  # ff %>% filter(KEY %in% (ff[duplicated(ff$KEY),]$KEY))
  
  # Deduplicate variants with same keys
  ff <- ff[(!duplicated(ff$KEY)),]
  ffpe <- ffpe[(!duplicated(ffpe$KEY)),]
  
  # Subset the variants by frequency if frequency argument is provided
  if (!is.null(var_freq)) {
    ff <- ff %>% filter(VAF > var_freq)
    ffpe <- ffpe %>% filter(VAF > var_freq)
  }
  
  # # Write the table with filtered variants
  # if (dim(ff)[1] != 0) {
  #   ff$SAMPLE_TYPE <- "FF"}
  # if (dim(ffpe)[1] != 0) {
  #   ffpe$SAMPLE_TYPE <- "FFPE"
  # }
  # all <- rbind(ff, ffpe)
  # all$class <- as.character(all$class)
  # write.table(all, file = paste0(ID, "_SNVs_", var_freq, ".tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
  
  # Remove RECURRENT variants (VAF < 0.1)
  ff <- ff %>% filter(!KEY %in% recurr_keys)
  ffpe <- ffpe %>% filter(!KEY %in% recurr_keys)
  
  # Summary table with concordance
  result <- data.frame(FF_TOTAL = dim(ff)[1],
                       FFPE_TOTAL = dim(ffpe)[1],
                       OVERLAP = sum(ff$KEY %in% ffpe$KEY), # bug
                       FF_UNIQ = ((dim(ff)[1]) - sum(ff$KEY %in% ffpe$KEY)),
                       FFPE_UNIQ = ((dim(ffpe)[1]) - sum(ff$KEY %in% ffpe$KEY)),
                       RECALL = ((sum(ff$KEY %in% ffpe$KEY))/(dim(ff)[1])),
                       PRECISION = (sum(ff$KEY %in% ffpe$KEY))/(dim(ffpe)[1]))
  
  return(result)
}

# Get all patient IDs
patientIDs <- as.character(unique(QC_portal_trios$PATIENT_ID))


# Run with different variant frequency thresholds
freq <- seq(from = 0, to = 0.3, by = 0.05)
lapply(freq, function(x){
  SNV_summary <- lapply(patientIDs, compareSNV_2, var_freq = x)
  SNV_summary <- bind_rows(SNV_summary)
  SNV_summary$PATIENT_ID <- patientIDs
  write.csv(SNV_summary, file = paste0("SNV_summary_62trios_noRecurr_", x, "_", today, ".csv"), row.names = F, quote = F)
})



