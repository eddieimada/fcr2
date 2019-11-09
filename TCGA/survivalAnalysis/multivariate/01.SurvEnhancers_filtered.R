###########################################################################
### Survival analysis by subtypes using e-lncRNA
### 
### Diego F Sanchez and Luigi Marchionni
###
###########################################################################
### First steps
### Set working directory
setwd(".")
#setwd("~/DropboxMech/FANTOM6/TCGA/survAnalysis/code/")

### Clean enviroments
rm(list=ls())

### Load libraries
library(Biobase)
require(ggplot2)
require(survival)
require(multtest)
require(edgeR)
#require(xtable)
#library(limma)
#require(edgeR)
#require(sva)
#library(dplyr)

###########################################################################
### Useful functions
###########################################################################

### Function that runs Cox analysis
runCoxTab <- function(dat, surv, pheno, tissue) {
  out <- t(apply(dat, 1, function(x, s, p, t) {
    if (t == "Breast") {
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage[stage == "stage x"] <- NA
      stage[stage %in% c("stage i", "stage ia", "stage ib", "stage ii", "stage iia", "stage iib")] <- 1
      stage[stage %in% c("stage iii", "stage iiia", "stage iiib", "stage iiic", "stage iv")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      ER <- as.character(pheno$xml_breast_carcinoma_estrogen_receptor_status)
      ER[ER=="Indeterminate"] <- NA
      PR <- as.character(pheno$xml_breast_carcinoma_progesterone_receptor_status)
      PR[PR=="Indeterminate"] <- NA
      HER <- as.character(pheno$xml_lab_proc_her2_neu_immunohistochemistry_receptor_status)
      HER[HER=="Indeterminate"] <- NA
      HER[HER=="Equivocal"] <- "Negative"
      cp <- coxph(surv ~ x + stage + age + ER + PR + HER)
    } else if (t == "Bile") {
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      CA19 <- pheno$xml_primary_pathology_ca_19_9_level
      gender <- factor(pheno$"cgc_case_gender")
      cp <- coxph(surv ~ x + stage + gender + age + CA19)
      
    } else if (t == "Bladder"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      # not reported      stage i     stage ii    stage iii     stage iv 
      #  2            4          134          149          144 
      stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      cp <- coxph(surv ~ x + stage + gender + age)
    } else if (t == "Colorectal"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage <- gsub("a$|b$|c$", "", stage)
      #stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      #stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      CEA <- pheno$xml_preoperative_pretreatment_cea_level
      # use CEA? 1/3 samples lost
      cp <- coxph(surv ~ x + stage + age + gender + CEA)
    } else if (t == "Esophagus"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage <- gsub("a$|b$|c$", "", stage)
      #stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      #stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      class <- pheno$xml_primary_pathology_histological_type
      # use CEA? 1/3 samples lost
      cp <- coxph(surv ~ x + stage + age + gender + class)
    } else if (t == "HeadNeck"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage <- gsub("a$|b$|c$", "", stage)
      #stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      #stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      PI <- pheno$xml_perineural_invasion_present
      smoke <- pheno$xml_number_pack_years_smoked
      # use CEA? 1/3 samples lost
      cp <- coxph(surv ~ x + stage + age + gender + PI + smoke)
    } else if (t == "Kidney"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage <- gsub("a$|b$|c$", "", stage)
      #stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      #stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      class <- pheno$gdc_cases.diagnoses.morphology
      cp <- coxph(surv ~ x + stage + age + gender + class)
    } else if (t == "Liver"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage <- gsub("a$|b$|c$", "", stage)
      #stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      #stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      cp <- coxph(surv ~ x + stage + gender + age)
    } else if (t == "Lung"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage <- gsub("a$|b$|c$", "", stage)
      #stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      #stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      class <- pheno$xml_icd_o_3_histology
      class <- as.numeric(gsub("/3", "", class))
      class[class == 8507] <- NA
      class <- ifelse(class < 8100, "SCC", "adenocarcinoma")
      smoke <- pheno$xml_number_pack_years_smoked
      cp <- coxph(surv ~ x + stage + age + gender + class + smoke)
    } else if (t == "Stomach"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage <- gsub("a$|b$|c$", "", stage)
      #stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      #stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      class <- pheno$xml_icd_o_3_histology
      class <- as.numeric(gsub("/3", "", class))
      class <- ifelse(class == 8145 | class == 8490, "Diffuse", "NOS")
      cp <- coxph(surv ~ x + stage + age + gender + class)
    } else if (t == "Thyroid"){
      gender <- factor(pheno$"cgc_case_gender")
      stage <- pheno$gdc_cases.diagnoses.tumor_stage
      stage[stage == "not reported"] <- NA
      stage <- gsub("a$|b$|c$", "", stage)
      #stage[stage %in% c("stage i", "stage ii", "stage iii")] <- 1
      #stage[stage %in% c("stage iv", "stage iva", "stage ivb")] <- 2
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      extension <- as.character(pheno$xml_extrathyroid_carcinoma_present_extension_status)
      extension<- ifelse(extension == "None", "None", "Ext")
      cp <- coxph(surv ~ x + stage + age + gender +  extension)
    } else if (t == "Uterus"){
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      #class <- pheno$gdc_cases.project.project_id
      grade <- as.character(pheno$xml_neoplasm_histologic_grade)
      grade <- ifelse(grade == "High Grade" | grade == "G3", 2, 1)
      grade[is.na(grade)] <- 1
      cp <- coxph(surv ~ x + age + grade)
    } else if (t == "Prostate"){
      age <- pheno$gdc_cases.diagnoses.days_to_birth
      age <- abs(age)/365
      gleason <- pheno$xml_stage_event_gleason_grading
      gleason <- substr(gleason,1,3)
      gleason <- as.numeric(substr(gleason,2,2))
      gleason[gleason == 0] <- 10
      gleason <- ifelse(gleason >= 4, 2, 1)
      cp <- coxph(surv ~ x + age + gleason)
    }
    else {
      cp <- coxph(surv ~ x )
    }
    test <- cox.zph(cp)
    test <- test$table[1,"p"]
    if(is.na(test)){
      out <- c("HR"=NA, 
               "low95ci"=NA,
               "up95ci"=NA,
               "Rsquare"=NA,
               "P-value"=NA)
    }
    else if (test <= 0.05) {
      out <- c("HR"=NA, 
               "low95ci"=NA,
               "up95ci"=NA,
               "Rsquare"=NA,
               "P-value"=NA)
    } else {
      cph <- summary(cp)
      out <- c("HR"=cph$conf.int[1,][1], 
               "low95ci"=cph$conf.int[1,][3],
               "up95ci"=cph$conf.int[1,][4],
               "Rsquare"=cph$rsq[[1]],
               "P-value"=coefficients(cph)["x", "Pr(>|z|)"])
    }}, s=surv, p=pheno,t=tissue))
  colnames(out) <- c("HR", "low95ci", "up95ci", "Rsquare", "P-value")
  myFDR <- mt.rawp2adjp(out[ , "P-value"], proc="BH")
  out <- cbind(out, myFDR$adjp[order(myFDR$index) , ])
  out <- out[ , -grep("rawp", colnames(out)) ]	
}


###########################################################################
### loading esets and making Cox analysis by types
###########################################################################
files <- list.files("../../DGE", recursive = TRUE,
                    pattern=".+elnc.+.rda$", full.names = TRUE)
### removing repeated objects or not desirable cancer types
files <- files[-grep("DGElist", files)]
# files <- files[-grep("Kidney/", files)]
# files <- files[-grep("Lung/", files)]
# files <- files[-grep("Uterus/", files)]
# files <- files[-grep("Colorectal/", files)]
files <- files[-grep("Rectum", files)]
files <- files[-grep("Colon/", files)]
files <- files[-grep("Kidney_sep/", files)]
files <- files[-grep("Lung_sep/", files)]
files <- files[-grep("Uterus_sep/", files)]

### allocating space for collecting objects
survObjectByTypes <- list()
survModelByTypes <- list()
coxModelByTypes <- list()

for(file in files) {

  load(file)
  type <- gsub(".+/", "", file)
  type <- gsub(".rda", "", type)
  organ <- gsub(".+/../DGE/", "", file)
  organ <- gsub("/.+", "", organ)
  eset <- get(type)
   
  ###########################################################################
  ### keeping primary cases and complete phenotypes
  
  ### extracting pheno data
  pheno <- pData(eset)
  ### keeping complete columns in pheno data
  keep.cols <- which(apply(pheno, 2, function(x) !all(is.na(x))))
  pheno <- pheno[,keep.cols]
  ### extracting expression data
  express <- exprs(eset)
  ### available sample types (normal, primary, mts, etc)
  table(pheno$gdc_cases.samples.sample_type)
  
  ### keeping primary tumor samples in Pheno data
  pheno <- pheno[which(pheno$gdc_cases.samples.sample_type == "Primary Tumor"), ]
  ## available sample types based on normal TCGA IDs
  table(pheno$gdc_cases.samples.sample_type)
 
   ## sample with frozen and ffpe could be available, delete ffpe
  idx <- which(pheno$gdc_cases.samples.is_ffpe == TRUE)
  if (length(idx) != 0){
    #express <- express[,-idx]
    pheno <- pheno[-idx,]
  }
  ## looking for cases with duplicated id after ffpe deleted
  dup <- unique(pheno$gdc_cases.submitter_id[duplicated(pheno$gdc_cases.submitter_id)])
  table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup] )
  
  
  ### Some cases have more than 1 primary tumor, DELETE ALL (NOT implemented)
  # tableNames <- names(table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup] ))
  # if (length(tableNames > 0)) {
  #   pheno <- pheno[-which(pheno$gdc_cases.submitter_id %in% tableNames), ]
  # }
  # table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup] )

  ### Some cases have more than 1 primary tumor, choose the first occurrence
  for (duplicates in dup) {
    x <- which(pheno$gdc_cases.submitter_id %in% duplicates)
    toDrop <- x[2:length(x)]
    pheno <- pheno[-toDrop, ]
  }
  table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup] )  
    
  ### Keeping primary tumor samples in Expression data
  sampleIDs <- rownames(pheno)
  express <- express[, colnames(express) %in% sampleIDs]
  ### Consistency check: Expression colnames are in the same order as pheno rownames
  all(rownames(pheno) == colnames(express))
  
  ### Removing eset object
  rm(list = type) #use list to remove the value of type and no type
  
  ###########################################################################
  ### Assembling survival object
  ## Combining time to event variables for surv_fu_days: 
  ## - Follow up time is in $gdc_cases.diagnoses.days_to_last_follow_up
  ## - Time to failure is not included in the previous (NA)
  ## - Time to failure for death is in $gdc_cases.diagnoses.days_to_death
  ## There are cgc_ and xml_ variables containing similar data with more NAs (less accurate?!)
  
  ### surv_fu_days captures the time to event
  pheno <- cbind(pheno, surv_fu_days = pheno$gdc_cases.diagnoses.days_to_last_follow_up)
  pheno$surv_fu_days[pheno$surv_fu_days == 0] <- NA # Changing uninformative time to event to "NA"
  ## merging days_to_last_follow_up with days_to_death
  which_dod <- which(!is.na(pheno$gdc_cases.diagnoses.days_to_death))
  pheno$surv_fu_days[which_dod] <- pheno$gdc_cases.diagnoses.days_to_death[which_dod]
  
  ### surv_deceased captures the event, any value other than NA represents dead
  pheno <- cbind(pheno, surv_deceased = pheno$gdc_cases.diagnoses.days_to_death)
  pheno$surv_deceased <- ifelse(is.na(pheno$surv_deceased), 0, 1) 
  
  ### Survival object
  survivalDOD <- Surv(pheno$surv_fu_days, event = pheno$surv_deceased)
  
  ### Survival model fitting: Kaplan-Meier estimator
  survModel <- survfit(survivalDOD ~ 1) 
  #summary(survModel)
  TCGAType <- gsub(".+-", "", pheno$gdc_cases.project.project_id[1])
  
  ### List object names using TCGA types or Organ types.
  #survObjectByTypes[[TCGAType]] <- survivalDOD
  survObjectByTypes[[organ]] <- survivalDOD
  #survModelByTypes[[TCGAType]] <- survModel
  survModelByTypes[[organ]] <- survModel
  
  ### Kaplan-Meier survival curve
  #plot(survModel)
  
  ### Filter low counts
  keepGns <- filterByExpr(express, group = pheno$surv_deceased, min.count = 5)
  table(keepGns)
  express <- log2(express[keepGns, ]+1)
  #following filter was not implemented. it is maintained for comparison.
  #keepGns2 <- rowSums(express > 5) >= 15
  #table(keepGns2)
  print(organ)
  ### Cox survival model fitting
  coxModel <- runCoxTab(dat=express, surv=survivalDOD, 
                        pheno=pheno, tissue=organ)
  
  ### List object names using TCGA types or Organ types.
  #coxModelByTypes[[TCGAType]] <- coxModel
  coxModelByTypes[[organ]] <- coxModel
}  

###########################################################################
### Saving object with Cox results by cancer type and session info
###########################################################################
save(coxModelByTypes, file="../objs/byOrganTCGATypesCoxelncRNA_filtered_log.rda")
save(survModelByTypes, file="../objs/byOrganTCGATypesSurvModels_filtered_log.rda")
save(survObjectByTypes, file="../objs/byOrganTCGATypesSurvObjects_filtered_log.rda")

### Session date and info
date()
sessionInfo()

### Clean
rm(list=ls())
quit(save = "no")
  
