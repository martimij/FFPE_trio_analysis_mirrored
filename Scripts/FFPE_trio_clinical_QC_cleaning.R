# FFPE trio analysis PART 1: Clinical and QC data cleaning
# Martina Mijuskovic
# May 2017



######### Load libraries  #########

library(dplyr)
library(ggplot2)
library(data.table)



######### Clean QC and clinical data for IIP trios #########

# Read QC data for all cancer samples ((Jan 26 2017))
QC_portal <- read.csv("./Data/ready_to_upload.01-26-17.csv")

# Summarize and check for dupliactes
dim(QC_portal)  # 1742 entries total
sum(duplicated(QC_portal))   # 0 duplicate rows
sum(duplicated(QC_portal$SAMPLE_WELL_ID))  # 0 duplicates
sum(duplicated(QC_portal$SAMPLE_LAB_ID))  # 81 duplicates

# Investigate SAMPLE_LAB_ID duplicates
QC_portal %>% filter(SAMPLE_LAB_ID %in% QC_portal$SAMPLE_LAB_ID[duplicated(QC_portal$SAMPLE_LAB_ID)])   # Duplicate entries for samples that didn't pass QC

# Any duplicate SAMPLE_LAB_IDs in samples that pass pre-sequencing QC? (NONE)
sum(duplicated(QC_portal %>% filter(SAMPLE_LAB_ID %in% QC_portal$SAMPLE_LAB_ID[duplicated(QC_portal$SAMPLE_LAB_ID)]) %>% filter(PASS_TO_SEQ == "PASS") %>% .$SAMPLE_LAB_ID))  # none

# Subset samples that pass QC, pass tumor and germline contamination (remove also those with no data available)
table(QC_portal$PASS_TO_SEQ, exclude = NULL)
# NOTE that I'm not filtering germline for contamination==PASS
QC_portal_passPreSeqQC <- QC_portal %>% filter(PASS_TO_SEQ %in% c("Pass", "TBD"), TUMOUR_CONTAMINATION != "Fail", GERMLINE_CONTAMINATION != "Fail")
table((QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE")) %>% .$TUMOUR_CONTAMINATION), exclude = NULL)
table((QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE == "GL") %>% .$GERMLINE_CONTAMINATION), exclude = NULL)
# Remove tumor samples (only) with no contamination data available
rmIDs <- QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE"), TUMOUR_CONTAMINATION == "") %>% .$SAMPLE_WELL_ID
QC_portal_passPreSeqQC <- QC_portal_passPreSeqQC %>% filter(!SAMPLE_WELL_ID %in% rmIDs)  # 1170 entries

# Flag trios
trios <- sapply(unique(QC_portal_passPreSeqQC$PATIENT_ID), function(x){
  if (sum(c("GL", "FFPE", "FF") %in% QC_portal_passPreSeqQC[QC_portal_passPreSeqQC$PATIENT_ID == x,]$SAMPLE_TYPE) == 3) {
    return(1)
  }
  else {
    return(0)
  }
})
names(trios) <- unique(QC_portal_passPreSeqQC$PATIENT_ID)

QC_portal_passPreSeqQC$Trio <- 0
QC_portal_passPreSeqQC[QC_portal_passPreSeqQC$PATIENT_ID %in% names(trios[trios == 1]),]$Trio <- 1

# Total number of trios that pass QC to sequencing
dim(QC_portal_passPreSeqQC %>% filter(Trio == 1))  # 86
table((QC_portal_passPreSeqQC %>% filter(Trio == 1, SAMPLE_TYPE %in% c("FF", "FFPE")) %>% .$TUMOUR_CONTAMINATION), exclude=NULL)

# Subset data for clean trios
QC_portal_trios <- QC_portal_passPreSeqQC %>% filter(Trio == 1)
table(table(QC_portal_trios$PATIENT_ID))  # NOTE that there are 3 samples for 28 patients but 2 patients have 4 samples
table(QC_portal_trios$SAMPLE_TYPE)  # There are 2 samples with extra germline sample
table((QC_portal_trios %>% filter(SAMPLE_TYPE =="GL") %>% .$GERMLINE_CONTAMINATION), exclude = NULL)
# Two patients with 2 GL entries, one entry in each doesn't Pass GERMLINE_CONTAMINATION
QC_portal_trios[QC_portal_trios$PATIENT_ID %in% (QC_portal_trios %>% filter(SAMPLE_TYPE == "GL", GERMLINE_CONTAMINATION != "Pass") %>% .$PATIENT_ID),]
# Remove extra GL entries
rmGLs <- QC_portal_trios %>% filter(SAMPLE_TYPE == "GL", GERMLINE_CONTAMINATION != "Pass") %>% .$SAMPLE_WELL_ID
QC_portal_trios <- QC_portal_trios %>% filter(!SAMPLE_WELL_ID %in% rmGLs)

### IMPORTANT: Remove any samples from WEST MIDLANDS, "RRK" (discontinued - unreliable)
QC_portal_trios <- QC_portal_trios %>% filter(CENTER_CODE != "RRK")

# Load IIP trios list from Alona to doublecheck
Alonas_trios <- read.table("./Data/trios.2016-01-26.txt", sep = "\t", colClasses = rep("character", 2))
names(Alonas_trios) <- c("PATIENT_ID", "STATUS")
Alonas_trios$ToUse <- 0
Alonas_trios[Alonas_trios$STATUS == "",]$ToUse <- 1
table(Alonas_trios$ToUse)  # 26 trios

# Check my trios against Alona's
trioIDs <- Alonas_trios %>% filter(ToUse == 1) %>% .$PATIENT_ID
sum(trioIDs %in% QC_portal_trios$PATIENT_ID)  # 26 (all there)

# Add tumor types
tumor_types <- read.table("./Data/cancer_disease_information_tum_2017-01-23_20-37-56.tsv", sep = "\t", header = T)
QC_portal_trios$TumorType <- tumor_types[match(QC_portal_trios$PATIENT_ID, tumor_types$participant_identifiers_id),]$disease_type_id
QC_portal_trios$TumorType <- as.character(QC_portal_trios$TumorType)

# Fix the unknown tumor type (from missing Oxford data)
missingTTid <- unique(QC_portal_trios %>% filter(TumorType == "Unknown") %>% .$PATIENT_ID)
missingOx <- read.csv("./Data/Oxford_samples_with_missing_data MC.csv", header = T)
QC_portal_trios[QC_portal_trios$PATIENT_ID == missingTTid,]$TumorType <- as.character(missingOx[missingOx$PATIENT_ID == missingTTid,]$X)

# # Read the current upload report, restrict to cancer, V4 and qc_passed
# today <- Sys.Date()
# system(paste0("wget ", "https://upload-reports.gel.zone/upload_report.", today, ".txt"))
# upload <- read.table(paste0("upload_report.", today, ".txt"), sep = "\t")
# colnames(upload) <- as.character(fread(paste0("upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
# upload <- upload %>% filter(`Delivery Version` == "V4", Status == "qc_passed", Type %in% c("cancer germline", "cancer tumour"))
# upload$Path <- as.character(upload$Path)

# Read an old upload report to get BAM paths
upload <- read.table("./Data/upload_report.2017-05-03.txt", sep = "\t")
colnames(upload) <- as.character(fread("./Data/upload_report.2017-05-03.txt", skip = 14, nrows = 1, header = F))
upload <- upload %>% filter(Platekey %in% QC_portal_trios$SAMPLE_WELL_ID, `Delivery Version` == "V4", Type %in% c("cancer germline", "cancer tumour"))
upload$Path <- as.character(upload$Path)

# Add BAM paths to the trios table
QC_portal_trios$BamPath <- upload[match(QC_portal_trios$SAMPLE_WELL_ID, upload$Platekey),]$Path

# # Manually add the BAM path for one sample (FFPE from 217000030, LP3000079-DNA_H01) - passes QC in v2 but not in v4, we are keeping it (NOT NECESSARY with above upload report, not missing)
# #upload %>% filter(Platekey == "LP3000079-DNA_H01") # Chose v4 path from unfiltered upload report
# QC_portal_trios[QC_portal_trios$SAMPLE_WELL_ID == "LP3000079-DNA_H01",]$BamPath <- "/genomes/by_date/2016-12-13/CANCP40747/CancerLP3000079-DNA_H01_NormalLP3000067-DNA_H08"

# Write clean trios QC + clinical data
write.csv(QC_portal_trios, file = "./Data/QC_IIP_FFPE_trios_clean.csv", quote = F, row.names = F)

# Cleanup
rm(Alonas_trios, missingOx, QC_portal, QC_portal_passPreSeqQC, tumor_types, missingTTid, rmGLs, rmIDs, trioIDs, trios)



######### Clean clinical data for BRC trios ######### 

# Read list of pilot samples (incl. BRC FFPE trios)
pilot <- read.csv("/Users/MartinaMijuskovic/FFPE_trio_analysis/Data/BRC_pilot_data/cancer_pilot_samples_David.csv")

# Read the two lists with tumor types
pilot_extra1 <- read.table("/Users/MartinaMijuskovic/FFPE_trio_analysis/Data/BRC_pilot_data/pilot_Matt_DB_data.txt", header = T)
pilot_extra2 <- read.table("/Users/MartinaMijuskovic/FFPE_trio_analysis/Data/BRC_pilot_data/Matt_DB_data.txt", header = T, sep = "\t")

# Explore
dim(pilot)  # 1014
table(pilot$study, exclude = NULL)
table(pilot$type, exclude = NULL)
table(pilot$type, exclude = NULL)
table(pilot$study, pilot$type)

# Reduce the list to BRC samples
pilot <- pilot %>% filter(study == "brc")
dim(pilot)  # 164

# Explore further
table(pilot$sex_problem, exclude = NULL)
table(pilot$modifier, pilot$type, exclude = NULL)

# Reduce to relevant columns
pilot <- pilot %>% dplyr::select(manifest_id, study, GelId, type, platekey, base_dir, modifier, assigned_diagnosis)

# Change names to fit 26 IIP trios
names(pilot)[3:5] <- c("PATIENT_ID", "SAMPLE_TYPE", "SAMPLE_WELL_ID")

# Flag trios
pilot$SAMPLE_TYPE <- as.character(pilot$SAMPLE_TYPE)
trios <- sapply(unique(pilot$PATIENT_ID), function(x){
  if (sum(c("GL", "FFPE", "FF") %in% pilot[pilot$PATIENT_ID == x,]$SAMPLE_TYPE) == 3) {
    return(1)
  }
  else {
    return(0)
  }
})
names(trios) <- unique(pilot$PATIENT_ID)
pilot$Trio <- 0
pilot[as.character(pilot$PATIENT_ID) %in% names(trios[trios == 1]),]$Trio <- 1

# Number of trios
sum(trios)  # 41
sum(pilot$Trio)/3  # 41.33333  (one FFPE sample extra!)
table(pilot$SAMPLE_TYPE, pilot$Trio, exclude = NULL)

# Find the origin of the extra sample
table(table(pilot$PATIENT_ID))  # 20 patients with 2 samples, 40 with trios, 1 with 4 samples
table(pilot$PATIENT_ID)  # 200000960 has 4 samples
sum(duplicated(pilot$SAMPLE_WELL_ID))  # 0
pilot %>% filter(PATIENT_ID == "200000960")  # this patient has 2 FFPE samples listed 
# (one may have failed but I have no info; older one?; LP2000830-DNA_E01 FFPE sample has only 20% purity)

# NOTE that one sample is denoted "EXPERIMENTAL" - best to remove; also missing diagnosis (PATIENT_ID==200000316)
pilot %>% filter(modifier == "EXPERIMENTAL")

# Explore extra data
table(pilot_extra1$CENTER, pilot_extra1$TUMOR, exclude = NULL)
table(pilot_extra1$CENTER, pilot_extra1$TYPE, exclude = NULL)
table(pilot_extra2$PILOT, pilot_extra2$TUMOR, exclude = NULL)
table(pilot_extra2$PILOT, pilot_extra2$type, exclude = NULL)

# Check if all trio samples have tumor type info
names(pilot_extra1)[1] <- "SAMPLE_WELL_ID"
names(pilot_extra2)[1] <- "SAMPLE_WELL_ID"
names(pilot_extra1)[5] <- "SAMPLE_TYPE"
names(pilot_extra2)[5] <- "SAMPLE_TYPE"
trio_ff_ids <- as.character(pilot %>% filter(Trio == 1, SAMPLE_TYPE == "FF") %>% .$SAMPLE_WELL_ID)
trio_ffpe_ids <- as.character(pilot %>% filter(Trio == 1, SAMPLE_TYPE == "FFPE") %>% .$SAMPLE_WELL_ID)
sum(!(trio_ff_ids %in% pilot_extra1[(!is.na(pilot_extra1$TUMOR)),]$SAMPLE_WELL_ID))  # 7 missing tumor info
sum(!(trio_ff_ids %in% pilot_extra2[(!is.na(pilot_extra2$TUMOR)),]$SAMPLE_WELL_ID))  # 5 missing tumor info
sum(!(trio_ffpe_ids %in% pilot_extra1[(!is.na(pilot_extra1$TUMOR)),]$SAMPLE_WELL_ID))  # 8 missing tumor info
sum(!(trio_ffpe_ids %in% pilot_extra2[(!is.na(pilot_extra2$TUMOR)),]$SAMPLE_WELL_ID))  # 6 missing tumor info

# Reduce extra tables to samples in the BRC pilot
pilot_extra1 <- pilot_extra1 %>% filter(SAMPLE_WELL_ID %in% pilot$SAMPLE_WELL_ID)  # this seems to be an incomplete subset of pilot_extra2
pilot_extra2 <- pilot_extra2 %>% filter(SAMPLE_WELL_ID %in% pilot$SAMPLE_WELL_ID)

# Check
sum(duplicated(pilot_extra2$SAMPLE_WELL_ID))  #0

# Check extra data for FF and FFPE trios
pilot_extra1 %>% filter(SAMPLE_WELL_ID %in% c(trio_ff_ids, trio_ffpe_ids))  # some marked "EXPT" (experimental?) under "PILOT"
pilot_extra2 %>% filter(SAMPLE_WELL_ID %in% c(trio_ff_ids, trio_ffpe_ids))  # some marked "EXPT" under "PILOT"

# Add tumor info from the 2nd file, clean up
pilot <- left_join(pilot, (pilot_extra2 %>% dplyr::select(SAMPLE_WELL_ID, TUMOR, SAMPLE_TYPE, CENTER, PILOT)))
table(pilot$Trio, pilot$TUMOR)
pilot[is.na(pilot$TUMOR),]$TUMOR <- "UNKNOWN"
pilot[pilot$TUMOR == "UNKOWN",]$TUMOR <- "UNKNOWN"
pilot$TUMOR <- as.character(pilot$TUMOR)

# PATIENT_IDs with missing tumor info
unique(pilot %>% filter(Trio == 1, TUMOR == "UNKNOWN") %>% .$PATIENT_ID)  # 200000926 200000929 200000953 200000954 200000960 200001299 200001310
pilot %>% filter(PATIENT_ID %in% (unique(pilot %>% filter(Trio == 1, TUMOR == "UNKNOWN") %>% .$PATIENT_ID)))

# Manually add tumor type using the dianogostic codes (http://www.icd10data.com/)
pilot[pilot$assigned_diagnosis %in% c("C50.9", "C50.3", "C50.8"),]$TUMOR <- "BREAST"
pilot[pilot$assigned_diagnosis == "C34.1",]$TUMOR <- "LUNG"
pilot[pilot$assigned_diagnosis %in% c("C18.7", "C19X", "C18.0"),]$TUMOR <- "COLORECTAL"

# Check tumor types
table(pilot$Trio, pilot$TUMOR, exclude = NULL) # No UNKNOWN tumor types in Trios left

# Subset to keep only trios
pilot <- pilot %>% filter(Trio == 1)  # 124

# Remove all "experimental" samples and the one that fails contamination check (200000336)
rm_ids <- as.character(c(unique(pilot %>% filter(PILOT == "EXPT") %>% .$PATIENT_ID), unique(pilot %>% filter(modifier == "EXPERIMENTAL") %>% .$PATIENT_ID), "200000336"))
pilot <- pilot %>% filter(!PATIENT_ID %in% rm_ids) # 108 (36 trios total)





######### Clean QC data for BRC trios ######### 

# Read QC data for BRC pilot samples
QC_BRC <- read.csv("/Users/MartinaMijuskovic/FFPE_trio_analysis/Data/BRC_pilot_data/BRC.QC.csv")

# Subset for BRC trios (NEW: 37 total)
QC_BRC$PATIENT_ID <- as.character(QC_BRC$PATIENT_ID)
QC_BRC <- QC_BRC %>% filter(PATIENT_ID %in% pilot$PATIENT_ID)
dim(QC_BRC)  # 111

# Clean names
names(QC_BRC)[2:3] <- c("SAMPLE_WELL_ID", "SAMPLE_TYPE")

# Merge pilot data with QC data
pilot$PATIENT_ID <- as.character(pilot$PATIENT_ID)
QC_BRC <- inner_join(pilot, QC_BRC)

# Change names to merge with other trios properly
names(QC_BRC)[6] <- "BamPath"
names(QC_BRC)[10] <- "TumorType"

# # Read the current upload report, restrict to cancer, V4 and qc_passed (connect to VPN!!!!)
# today <- Sys.Date()
# system(paste0("wget ", "https://upload-reports.gel.zone/upload_report.", today, ".txt"))
# upload <- read.table(paste0("upload_report.", today, ".txt"), sep = "\t")
# colnames(upload) <- as.character(fread(paste0("upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
# upload <- upload %>% filter(Platekey %in% QC_BRC$SAMPLE_WELL_ID, `Delivery Version` == "V4", Type %in% c("cancer germline", "cancer tumour"))
# upload$Path <- as.character(upload$Path)

# Read the BAM paths from an upload report
upload <- read.table("./Data/upload_report.2017-05-03.txt", sep = "\t")
colnames(upload) <- as.character(fread("./Data/upload_report.2017-05-03.txt", skip = 14, nrows = 1, header = F))
upload <- upload %>% filter(Platekey %in% QC_BRC$SAMPLE_WELL_ID, `Delivery Version` == "V4", Type %in% c("cancer germline", "cancer tumour"))
upload$Path <- as.character(upload$Path)

# Add new (v4) BAM paths to the BRC table
QC_BRC$BamPath <- upload[match(QC_BRC$SAMPLE_WELL_ID, upload$Platekey),]$Path

# List BRC trio samples with non-pass status
upload %>% filter(Platekey %in% QC_BRC$SAMPLE_WELL_ID, Status != "qc_passed") %>% dplyr::select(Platekey, Type, DeliveryID, `Delivery Date`, `Delivery Version`, `BAM Date`, Status) # none

# Remove empty columns
table(QC_BRC$TUMOUR_TYPE, QC_BRC$COLLECTING_DATE, exclude = NULL)
QC_BRC <- QC_BRC %>% select(-(TUMOUR_TYPE), -(COLLECTING_DATE))

# Write clean BRC QC + clinical data
write.csv(QC_BRC, file = "./Data/QC_BRC_FFPE_trios_clean.csv", quote = F, row.names = F)

# Cleanup
rm(pilot, pilot_extra1, pilot_extra2, rm_ids, trio_ff_ids, trio_ffpe_ids, trios)



######### Merge BRC and IIP trios #########


# # Read QC data for 26 IIP trios (if not loaded)
# QC_portal_trios <- read.csv("./Data/QC_IIP_FFPE_trios_clean.csv")
# 
# # Read QC data for 37 BRC trios (if not loaded)
# QC_BRC <- read.csv("./Data/QC_BRC_FFPE_trios_clean.csv")

# Change tumor type entries to all CAPS
QC_portal_trios$TumorType <- toupper(QC_portal_trios$TumorType)

# Merge BRC with initial 26 trios
QC_portal_trios$PATIENT_ID <- as.character(QC_portal_trios$PATIENT_ID)
QC_portal_trios <- bind_rows(QC_portal_trios, QC_BRC)

# Sanity check
dim(QC_portal_trios)[1]/3  # 62 trios (as expected)
sum(duplicated(QC_portal_trios$SAMPLE_WELL_ID))  # 0
summary(QC_portal_trios)  # tumor purity missing for BRC samples


# Write out the full FFPE trio QC data table
write.csv(QC_portal_trios, file = "./Data/QC_all62_FFPE_trios_clean.csv", quote = F, row.names = F)







