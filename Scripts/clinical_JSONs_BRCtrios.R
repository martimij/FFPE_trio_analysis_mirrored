# Create JSON files with clinical information for BRC trios
# Martina Mijuskovic
# May 2017

############  Create json files with clinical info for BRC trios ############  

# Read example json with clinical info
clin_json <- jsonlite::fromJSON("/Users/MartinaMijuskovic/FFPE/217000058_clinical_info_example.json")
#clin_json <- fromJSON(file = "217000058_clinical_info_example.json")

# Modify and write with patient ID in the file name
clinJsonModify <- function(patient_ID){
  
  # Read example json with clinical info
  clin_json <- jsonlite::fromJSON("/Users/MartinaMijuskovic/FFPE/217000058_clinical_info_example.json")
  
  ### Modify cancerSamples table (tumorType=TumorType, labId=NA, sampleType=(tumor, germline), 
  # sampleId=SAMPLE_WELL_ID, gelPhase=OXFORD, tumorSubType=NA, tumorContent=NA, preservationMethod=SAMPLE_TYPE)
  clin_json$cancerSamples$sampleId <- BRC_trios[BRC_trios$PATIENT_ID == patient_ID,]$SAMPLE_WELL_ID
  clin_json$cancerSamples$labId <- NA
  clin_json$cancerSamples$tumorSubType <- NA
  clin_json$cancerSamples$tumorContent <- NA
  clin_json$cancerSamples$gelPhase <- "OXFORD"
  clin_json$cancerSamples$tumorType <- BRC_trios[match(clin_json$cancerSamples$sampleId, BRC_trios$SAMPLE_WELL_ID),]$TumorType
  clin_json$cancerSamples$preservationMethod <- BRC_trios[match(clin_json$cancerSamples$sampleId, BRC_trios$SAMPLE_WELL_ID),]$SAMPLE_TYPE
  clin_json$cancerSamples[clin_json$cancerSamples$preservationMethod == "GL",]$tumorType <- NA
  clin_json$cancerSamples[clin_json$cancerSamples$preservationMethod == "GL",]$sampleType <- "germline"
  clin_json$cancerSamples[clin_json$cancerSamples$preservationMethod %in% c("FF", "FFPE"),]$sampleType <- "tumor"
  
  ### Modify cancerDemographics list
  #clin_json$cancerDemographics <- as.data.frame(clin_json$cancerDemographics)
  clin_json$cancerDemographics$center <- "OXFORD"
  clin_json$cancerDemographics$labkeyParticipantId <- patient_ID
  clin_json$cancerDemographics$sex <- NA
  clin_json$cancerDemographics$centerPatientId <- patient_ID
  clin_json$cancerDemographics$gelId <- patient_ID
  
  ### Modify matchedSamples
  clin_json$matchedSamples$germlineSampleId <- as.character(rep((clin_json$cancerSamples[clin_json$cancerSamples$preservationMethod == "GL",]$sampleId), 2))
  clin_json$matchedSamples$tumorSampleId <- as.character(clin_json$cancerSamples[clin_json$cancerSamples$sampleType == "tumor",]$sampleId)
  
  ### Change NA formats (NAs will be converted back to "null", but if it stays "null" it will be converted to an empty dictionary)
  clin_json$cancerDemographics$sampleId <- NA
  clin_json$cancerDemographics$additionalInformation <- NA
  clin_json$cancerDemographics$primaryDiagnosis <- NA
  clin_json$cancerDemographics$assignedICD10 <- NA
  
  ### Write json (!!!! use auto_unbox=T to prevent single element chr vectors to end up as lists in json)
  write(jsonlite::toJSON(clin_json, na = "null", pretty = T, auto_unbox = T), file = paste0(patient_ID, "_clin.json"))
}

patientIDs <- as.character(unique(BRC_trios$PATIENT_ID))

# Change directory
setwd("/Users/MartinaMijuskovic/FFPE/BRC_jsons")

# Create clinical jsons for all 36 BRC trios
sapply(patientIDs, clinJsonModify)


############  Add one extra BRC sample (failed QC) ############  

# Adding two samples that failed ContEst: LP2000514-DNA_A01 and LP2000499-DNA_A01 
# We will process them to look at their behavior

### Re-creating QC table with BRC trios

# Read list of pilot samples (incl. BRC FFPE trios)
pilot <- read.csv("/Users/MartinaMijuskovic/Documents/FFPE/Pilot data/cancer_pilot_samples_David.csv")

# Read the two lists with tumor types
pilot_extra1 <- read.table("/Users/MartinaMijuskovic/Documents/FFPE/Pilot data/pilot_Matt_DB_data.txt", header = T)
pilot_extra2 <- read.table("/Users/MartinaMijuskovic/Documents/FFPE/Pilot data/Matt_DB_data.txt", header = T, sep = "\t")

# Reduce the list to BRC samples
pilot <- pilot %>% filter(study == "brc")
dim(pilot)  # 164

# Reduce
pilot <- pilot %>% dplyr::select(manifest_id, study, GelId, type, platekey, base_dir, modifier, assigned_diagnosis)

# Change names to fit 26 previous trios analysis
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
sum(pilot$Trio)/3  # 41.33333  (one FFPE sample extra!) # 200000960 has 4 samples
table(pilot$SAMPLE_TYPE, pilot$Trio, exclude = NULL)

# NOTE that one sample is denoted "EXPERIMENTAL" - best to remove; also missing diagnosis (PATIENT_ID==200000316)
pilot %>% filter(modifier == "EXPERIMENTAL")

# Check if all trio samples have tumor type info
names(pilot_extra1)[1] <- "SAMPLE_WELL_ID"
names(pilot_extra2)[1] <- "SAMPLE_WELL_ID"
names(pilot_extra1)[5] <- "SAMPLE_TYPE"
names(pilot_extra2)[5] <- "SAMPLE_TYPE"
trio_ff_ids <- as.character(pilot %>% filter(Trio == 1, SAMPLE_TYPE == "FF") %>% .$SAMPLE_WELL_ID)
trio_ffpe_ids <- as.character(pilot %>% filter(Trio == 1, SAMPLE_TYPE == "FFPE") %>% .$SAMPLE_WELL_ID)

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

# NEW: Before removing any samples, look for two that fail ContEst
pilot %>% filter(SAMPLE_WELL_ID %in% c("LP2000514-DNA_A01", "LP2000499-DNA_A01"))  # Both of these are from same patient, FF and FFPE (200000336)

# Remove all "experimental" samples BUT KEEP patient that fails contamination check (200000336)
rm_ids <- as.character(c(unique(pilot %>% filter(PILOT == "EXPT") %>% .$PATIENT_ID), unique(pilot %>% filter(modifier == "EXPERIMENTAL") %>% .$PATIENT_ID)))
pilot <- pilot %>% filter(!PATIENT_ID %in% rm_ids)

# Read QC data for BRC pilot samples
QC_BRC <- read.csv("/Users/MartinaMijuskovic/Documents/FFPE/Pilot data/BRC.QC.csv")

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

# Read the current upload report, restrict to cancer, V4 and qc_passed (connect to VPN!!!!)
today <- Sys.Date()
system(paste0("wget ", "https://upload-reports.gel.zone/upload_report.", today, ".txt"))
upload <- read.table(paste0("upload_report.", today, ".txt"), sep = "\t")
colnames(upload) <- as.character(fread(paste0("upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
upload <- upload %>% filter(Platekey %in% QC_BRC$SAMPLE_WELL_ID, `Delivery Version` == "V4", Type %in% c("cancer germline", "cancer tumour"))
upload$Path <- as.character(upload$Path)

# Add new (v4) BAM paths to the BRC table
QC_BRC$BamPath <- upload[match(QC_BRC$SAMPLE_WELL_ID, upload$Platekey),]$Path

# List BRC trio samples with non-pass status
upload %>% filter(Platekey %in% QC_BRC$SAMPLE_WELL_ID, Status != "qc_passed") %>% dplyr::select(Platekey, Type, DeliveryID, `Delivery Date`, `Delivery Version`, `BAM Date`, Status)
QC_BRC %>% filter(SAMPLE_WELL_ID == "LP2000485-DNA_A01")  # germline of patient 200000336


### Add BRC data to initial 26 trios

# Read QC data of initial 26 trios
QC_portal_trios <- read.csv("QC_portal_trios_final.csv")

# Change tumor type entries to all CAPS
QC_portal_trios$TumorType <- toupper(QC_portal_trios$TumorType)

# Merge BRC with initial 26 trios
QC_portal_trios$PATIENT_ID <- as.character(QC_portal_trios$PATIENT_ID)
QC_portal_trios <- bind_rows(QC_portal_trios, QC_BRC)

# Write out the full FFPE trio QC data table
write.csv(QC_portal_trios, file = "QC_portal_63_trios.csv", quote = F, row.names = F)  # NOTE that no QC info is available for contaminated tumor samples from 200000336

# Write the clinical json for 200000336
#QC_portal_trios <- read.csv("/Users/MartinaMijuskovic/FFPE/QC_portal_63_trios.csv")
BRC_trios <- QC_portal_trios %>% filter(study == "brc")
setwd("/Users/MartinaMijuskovic/FFPE/BRC_jsons")
clinJsonModify("200000336")


############# Add invented LabID to json files ############

# The invented labID will be 10 digit number, consisting of 9-digit patientID and 1,2,3

QC_portal_trios <- read.csv("/Users/MartinaMijuskovic/FFPE/QC_portal_63_trios.csv")
BRC_trios <- QC_portal_trios %>% filter(study == "brc")

# Modify the function that creates clinical json files to include LabID
clinJsonModify2 <- function(patient_ID){
  
  # Read example json with clinical info
  clin_json <- jsonlite::fromJSON("/Users/MartinaMijuskovic/FFPE/217000058_clinical_info_example.json")
  
  ### Modify cancerSamples table (tumorType=TumorType, labId=10digitN, sampleType=(tumor, germline), 
  # sampleId=SAMPLE_WELL_ID, gelPhase=OXFORD, tumorSubType=NA, tumorContent=NA, preservationMethod=SAMPLE_TYPE)
  clin_json$cancerSamples$sampleId <- BRC_trios[BRC_trios$PATIENT_ID == patient_ID,]$SAMPLE_WELL_ID
  clin_json$cancerSamples$labId <- paste0(patient_ID, 1:3)
  clin_json$cancerSamples$tumorSubType <- NA
  clin_json$cancerSamples$tumorContent <- NA
  clin_json$cancerSamples$gelPhase <- "OXFORD"
  clin_json$cancerSamples$tumorType <- BRC_trios[match(clin_json$cancerSamples$sampleId, BRC_trios$SAMPLE_WELL_ID),]$TumorType
  clin_json$cancerSamples$preservationMethod <- BRC_trios[match(clin_json$cancerSamples$sampleId, BRC_trios$SAMPLE_WELL_ID),]$SAMPLE_TYPE
  clin_json$cancerSamples[clin_json$cancerSamples$preservationMethod == "GL",]$tumorType <- NA
  clin_json$cancerSamples[clin_json$cancerSamples$preservationMethod == "GL",]$sampleType <- "germline"
  clin_json$cancerSamples[clin_json$cancerSamples$preservationMethod %in% c("FF", "FFPE"),]$sampleType <- "tumor"
  
  ### Modify cancerDemographics list
  #clin_json$cancerDemographics <- as.data.frame(clin_json$cancerDemographics)
  clin_json$cancerDemographics$center <- "OXFORD"
  clin_json$cancerDemographics$labkeyParticipantId <- patient_ID
  clin_json$cancerDemographics$sex <- NA
  clin_json$cancerDemographics$centerPatientId <- patient_ID
  clin_json$cancerDemographics$gelId <- patient_ID
  
  ### Modify matchedSamples
  clin_json$matchedSamples$germlineSampleId <- as.character(rep((clin_json$cancerSamples[clin_json$cancerSamples$preservationMethod == "GL",]$sampleId), 2))
  clin_json$matchedSamples$tumorSampleId <- as.character(clin_json$cancerSamples[clin_json$cancerSamples$sampleType == "tumor",]$sampleId)
  
  ### Change NA formats (NAs will be converted back to "null", but if it stays "null" it will be converted to an empty dictionary)
  clin_json$cancerDemographics$sampleId <- NA
  clin_json$cancerDemographics$additionalInformation <- NA
  clin_json$cancerDemographics$primaryDiagnosis <- NA
  clin_json$cancerDemographics$assignedICD10 <- NA
  
  ### Write json (!!!! use auto_unbox=T to prevent single element chr vectors to end up as lists in json)
  write(jsonlite::toJSON(clin_json, na = "null", pretty = T, auto_unbox = T), file = paste0(patient_ID, "_clin.json"))
}

# Re-run to create clinical json files for BRC samples, including the LabID
patientIDs <- as.character(unique(BRC_trios$PATIENT_ID))

# Change directory
setwd("/Users/MartinaMijuskovic/FFPE/BRC_jsons")

# Create clinical jsons for all 36 BRC trios
sapply(patientIDs, clinJsonModify2)


