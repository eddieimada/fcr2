###########################################################################
### Survival analysis by subtypes using e-lncRNA
### 
### Diego F Sanchez and Luigi Marchionni
###
###########################################################################
### First steps
### Set working directory
setwd(".")
#/Users/diegosanchez/DropboxMech/FANTOM6/TCGA/survAnalysis/code
### Clean enviroments
rm(list=ls())

### Load libraries
library(Biobase)
require(ggplot2)
require(survival)
require(multtest)
#require(xtable)
#library(limma)
#require(edgeR)
#require(sva)
#library(dplyr)

###########################################################################
### Useful functions
###########################################################################

### Function that runs Cox analysis
runCoxTab <- function(dat, surv) {
  out <- t(apply(dat, 1, function(x, s) {
    cph <- summary(coxph(surv ~ x))
    out <- c("HR"=cph$conf.int[1], 
             "low95ci"=cph$conf.int[3],
             "up95ci"=cph$conf.int[4],
             "Rsquare"=cph$rsq[[1]],
             "P-value"=coefficients(cph)["x", "Pr(>|z|)"])
  }, s=surv))
  myFDR <- mt.rawp2adjp(out[ , "P-value"], proc="BH")
  out <- cbind(out, myFDR$adjp[order(myFDR$index) , ])
  out <- out[ , -grep("rawp", colnames(out)) ]	
}


###########################################################################
### loading esets and making Cox analysis by types
###########################################################################
files <- list.files("../../DGE/", recursive = TRUE,
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
  organ <- gsub(".+//", "", file)
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
  
  ### Cox survival model fitting
  coxModel <- runCoxTab(dat=express, surv=survivalDOD)
  
  ### List object names using TCGA types or Organ types.
  #coxModelByTypes[[TCGAType]] <- coxModel
  coxModelByTypes[[organ]] <- coxModel
}  

###########################################################################
### Saving object with Cox results by cancer type and session info
###########################################################################
save(coxModelByTypes, file="../objs/byOrganTCGATypesCoxelncRNA.rda")
save(survModelByTypes, file="../objs/byOrganTCGATypesSurvModels.rda")
save(survObjectByTypes, file="../objs/byOrganTCGATypesSurvObjects.rda")

### Session date and info
date()
sessionInfo()

### Clean
rm(list=ls())
quit(save = "no")
  