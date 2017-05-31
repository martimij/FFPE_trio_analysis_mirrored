# Create JSON files with clinical information for BRC trios
# Martina Mijuskovic
# May 2017

library(jsonlite)


# Read in the QC data for trios, subset to BRC trios
QC_portal_trios <- read.csv("/Users/MartinaMijuskovic/FFPE_trio_analysis/Data/QC_all62_FFPE_trios_clean.csv")
BRC_trios <- QC_portal_trios %>% filter(study == "brc")

# Modify the function that creates clinical json files to include LabID
clinJsonModify <- function(patient_ID){
  
  # Read example json with clinical info
  clin_json <- jsonlite::fromJSON("/Users/MartinaMijuskovic/FFPE/217000058_clinical_info_example.json")
  
  ### Modify cancerSamples table (tumorType=TumorType, labId=10digitN, sampleType=(tumor, germline), 
  # sampleId=SAMPLE_WELL_ID, gelPhase=OXFORD, tumorSubType=NA, tumorContent=NA, preservationMethod=SAMPLE_TYPE)
  clin_json$cancerSamples$sampleId <- BRC_trios[BRC_trios$PATIENT_ID == patient_ID,]$SAMPLE_WELL_ID
  clin_json$cancerSamples$labId <- paste0(patient_ID, 1:3)  # invented labId (9-digit patientID and 1,2,3)
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
setwd("./Data/BRC_jsons/")

# Create clinical jsons for all 36 BRC trios
sapply(patientIDs, clinJsonModify)


