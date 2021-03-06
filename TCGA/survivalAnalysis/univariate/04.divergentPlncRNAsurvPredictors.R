###########################################################################
### Survival analysis by subtypes using divergent-plncRNA
### 
### Diego F Sanchez
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

### Function to save data.frames from a list to files 
saveDFtoFiles <- function(x, y) {
  # x = names(list), y = list
  temp.df <- y[[x]]
  temp.df <- temp.df[ order(abs(temp.df$HR), decreasing=TRUE), ]
  temp.df <- temp.df[ order(temp.df$BH), ]
  write.csv(temp.df, file = paste("../text/divergentPlncRNA/SigDivergentPlncRNA_", x, ".csv", sep = ""))
}

###########################################################################
### loading esets and making Cox analysis by types
###########################################################################
files <- list.files("../../DGE/", recursive = TRUE,
                    pattern=".+divergent.+.rda$", full.names = TRUE)
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
  
  ########## WARNING ###########
  ### Some data do not fit using coxph. Expression level for these problematic 
  # points are removed here:
  # For Breast gene #4493 generates NA/NaN/Inf in foreign function call (arg 6) error (coxph)
  if(grepl("Breast", file)){
    # remove expression for gene #4493
    express <- express[-4493,]
  }
  ##############################
  
  ### Cox survival model fitting
  coxModel <- runCoxTab(dat=express, surv=survivalDOD)
  
  ### List object names using TCGA types or Organ types.
  #coxModelByTypes[[TCGAType]] <- coxModel
  coxModelByTypes[[organ]] <- coxModel
}  

###########################################################################
### Number of cases with survival time and adverse events
survivalTable <- t(sapply(survModelByTypes, function (x) {
  temp <- summary(x)$table
  survivalTable <- c(temp[2], temp[4], temp[7])
}))

### Check table
class(survivalTable)
survivalTable

### renaming the matrix
colnames(survivalTable) <- c("Cases","Events","Median time") #renaming the matrix
survivalTable

###########################################################################
### Number of significant e-lncRNA by univariate Cox
countBySignificance <- t(sapply(coxModelByTypes, function (x) {table(x[, "BH"] <= 0.05)}))
### Check if there are some cancer types with no predictors
class(countBySignificance)
typeof(countBySignificance)
countBySignificance

### It is not a numeric matrix because some columns were missed by not T or F values
#making vectors
tempSignificant <- lapply(lapply(countBySignificance, unlist), "[", c("FALSE", "TRUE")) 
#making the numeric matrix (merging previous values by columns)
countBySignificance <- do.call(rbind, tempSignificant) 
#renaming the matrix rows and columns
colnames(countBySignificance) <- c("Non-significant", "FDR <= 0.05") 
rownames(countBySignificance) <- names(coxModelByTypes)
#change NAs to "0"
countBySignificance[ is.na(countBySignificance)] <- 0
#checking the final matrix
countBySignificance

### Making a data.frame with available data from models
# binding previous matrices
summaryDF <- as.data.frame(cbind(countBySignificance, survivalTable))
# sorting by FDR
summaryDF <- summaryDF[order(summaryDF$`FDR <= 0.05`, decreasing = TRUE), ]
# creating types and renaming rows
Types <- factor(row.names(summaryDF), levels = row.names(summaryDF))
summaryDF <- cbind(Types, summaryDF)
row.names(summaryDF) <- 1:nrow(summaryDF)
# checking final data.frame
summaryDF
# Saving into significantDivergentPromByTCGAtypes.csv
write.csv(summaryDF, file = "../text/divergentPlncRNA/significantDivergentPromByTCGAtypes.csv")

###########################################################################
### data.frames with Significant genes (BH <= 0.05)
significantDivergentProm <- lapply(coxModelByTypes, function(x){
  tmp <- as.data.frame(x[!is.na(x[, "HR"]), ]) # It is generating NAs because there are NAs in the Cox Models.
  tmp <- tmp[tmp[,"BH"] <= 0.05, ]
})

### Saving data.frames with elncRNAs by TCGA types
invisible(lapply(names(significantDivergentProm), saveDFtoFiles, y = significantDivergentProm))

##########################################################################
### Common genes among subtypes
# Retaining just names
dpromNames <- lapply(significantDivergentProm, row.names)
dpromNames <- dpromNames[sapply(dpromNames, length) > 0]

# elncRNAs shared among all cancer types
dpromNamesCommon <- Reduce(f="intersect", x = dpromNames)
length(dpromNamesCommon)
dpromNamesCommon

# Predictive elncRNAs
dpromNamesPredictives <- Reduce(f="union", x = dpromNames)
length(dpromNamesPredictives)
write.csv(x = dpromNamesPredictives, file = "../text/divergentPlncRNA/union_of_predictive_dplncRNAs_through_tumors_list.csv")
#dpromNamesPredictives

###########################################################################
### closing and saving session info
###########################################################################

### Session date and info
date()
sessionInfo()

### Clean
rm(list=ls())
quit(save = "no")


