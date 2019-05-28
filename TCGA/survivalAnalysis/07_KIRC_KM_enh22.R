###########################################################################
### Survival analysis by TCGA cancer types using e-lncRNA
### KIRC survival curve for Enhancer 22
### Diego F Sanchez and Luigi Marchionni
###
###########################################################################
### Set working directory
setwd(".")

### Clean enviroments
rm(list=ls())

### Load libraries
library(Biobase)
#require(VennDiagram)
require(ggplot2)
require(survival)
#require(xtable)
require(multtest)
#require(plyr)
require(gtools)
require(survminer)

### ONLY KIDNEY DATA LOOK PROMISING. 
### Preparing Kidney cancer data
### Loading Kidney data
load("~/DropboxMech/FANTOM6/TCGA/DGE/Kidney_sep/objs/KIRCelncRNA.rda")
### extracting pheno data
pheno <- pData(KIRCelncRNA)
### keeping complete columns in pheno data
keep.cols <- which(apply(pheno, 2, function(x) !all(is.na(x))))
pheno <- pheno[,keep.cols]
### extracting expression data
express <- exprs(KIRCelncRNA)
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


### Creating prognostic variable for ENSG00000272666
ENSG00000272666values <- express[rownames(express) == "ENSG00000272666", ]
ENSG00000272666 <- ifelse(ENSG00000272666values > median(ENSG00000272666values), 1, 0)

### Selecting Kidney survival object
KidneySurv <- survivalDOD

### Fitting survival model according ENSG00000272666
kidneyEnhancerSurv <- survfit(KidneySurv ~ ENSG00000272666)

### Drawing Kaplan-Meier curve
jpeg(filename = "../figs/KIRC_KM_curve.jpg", width = 4000, height = 4000, res = 300)
ggsurvkm <- ggsurvplot(kidneyEnhancerSurv, data = KidneySurv, pval = TRUE, pval.method = TRUE,
                       xlab = "Overall survival in days", ylab = "Survival probability", pval.size = 10,
                       legend.lab = c("Low-risk, n = 266", "High-risk, n = 265"), legend = c(0.75,0.8),
                       palette = c("blue", "red"), size = 2,
                       #font.x = 5, font.y = 5, 
                       pval.coord = c(1200 ,0.1), pval.method.coord = c(500, 0.1))
ggsurvkm$plot <- ggsurvkm$plot + 
  theme(
    axis.text.x=element_text(color = "black", size=30),
    axis.text.y=element_text(color = "black", size=30),
    axis.title.x=element_text(size = 40),
    axis.title.y=element_text(size = 40),
    legend.text=element_text(size = 34),
    legend.title=element_text(size = 34),
    legend.key.size = unit(1.5, "cm"))
ggsurvkm
dev.off()
