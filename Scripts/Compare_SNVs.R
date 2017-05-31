# Martina Mijuskovic
# FFPE project
# Compares SNV calls between FFPE and FF samples from FFPE trios (HPC)

# NOTE: Installing R packages on HPC: use lib = "~/R/x86_64-pc-linux-gnu-library/3.3"

library(dplyr)
library(jsonlite)

# Working directory on the HPC
setwd("/home/mmijuskovic/FFPE_trio_analysis/SNV_trio_comparison")

# Today's date
today <- Sys.Date()



############# Get SNV VCF and JSON paths ############# 

# Load the manifest (HPC)
QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE_trio_analysis/QC_all62_FFPE_trios_clean.csv")  # all 62 trios
# QC_portal_trios <- read.csv("./Data/QC_portal_trios_final.csv") # local

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
  
  # Get json paths
  ff_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FF") %>% .$json_path
  ffpe_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FFPE") %>% .$json_path
  
  ### Read json files: get chromosome, position, VAF, tier, etc
  # FF
  #ff_json <- fromJSON("LP2000907-DNA_F02_tiering.json", flatten = T) 
  ff_json <- fromJSON(ff_path, flatten = T)
  ff <- ff_json %>% dplyr::select(
    reportedVariantCancer.chromosome, reportedVariantCancer.position, reportedVariantCancer.reference, reportedVariantCancer.alternate, 
    reportedVariantCancer.VAF)
  ff$tier <- sapply(1:dim(ff)[1], function(x){
    ff_json$reportedVariantCancer.reportEvents[[x]]$tier})
  ff$class <- sapply(1:dim(ff)[1], function(x){
    ff_json$reportedVariantCancer.reportEvents[[x]]$soNames})
  names(ff) <- c("chr", "pos", "ref", "alt", "VAF", "tier", "class")
  ff$KEY <- sapply(1:dim(ff)[1], function(x){
    paste(ff$chr[x], ff$pos[x], ff$ref[x], ff$alt[x], sep = "_")
  })

  # FFPE
  #ffpe_json <- fromJSON("LP2000906-DNA_H01_tiering.json", flatten = T)
  ffpe_json <- fromJSON(ffpe_path, flatten = T)
  ffpe <- ffpe_json %>% dplyr::select(
    reportedVariantCancer.chromosome, reportedVariantCancer.position, reportedVariantCancer.reference, reportedVariantCancer.alternate, 
    reportedVariantCancer.VAF)
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
  
  # Summary table with concordance
  result <- data.frame(FF_TOTAL = dim(ff)[1],
                       FFPE_TOTAL = dim(ffpe)[1],
                       OVERLAP = sum(ff$KEY %in% ffpe$KEY), # bug
                       FF_UNIQ = ((dim(ff)[1]) - sum(ff$KEY %in% ffpe$KEY)),
                       FFPE_UNIQ = ((dim(ffpe)[1]) - sum(ff$KEY %in% ffpe$KEY)),
                       RECALL = ((sum(ff$KEY %in% ffpe$KEY))/(dim(ff)[1])),
                       PRECISION = (sum(ff$KEY %in% ffpe$KEY))/(dim(ffpe)[1]))

 #write.table(result, file = paste0(patientID, "_SNV_concord", ".tsv"), row.names = F, quote = F, sep = "\t")
 return(result)
}

# Get all patient IDs
patientIDs <- as.character(unique(QC_portal_trios$PATIENT_ID))

# Run SNV comparison for each patient ID (all allele frequencies)
SNV_summary <- lapply(patientIDs, compareSNV)  # errors :/ " no applicable method for 'select_' applied to an object of class "list""
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


