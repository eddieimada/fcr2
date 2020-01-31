###########################################################################
### Survival analysis by TCGA cancer types
### Kidney survival curves by RNA types
### Diego F Sanchez
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
#require(ggpubr)

###########################################################################
### Preparing Kidney cancer data
### Loading Kidney data
load("~/DropboxMech/FANTOM6/TCGA/DGE/Kidney/objs/KidneyCodingmRNA.rda")
### extracting pheno data
pheno <- pData(KidneyCodingmRNA)
### keeping complete columns in pheno data
keep.cols <- which(apply(pheno, 2, function(x) !all(is.na(x))))
pheno <- pheno[,keep.cols]
### extracting expression data
express <- exprs(KidneyCodingmRNA)
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

###########################################################################
### Creating prognostic variable for ENSG00000065328
ENSG00000065328values <- express[rownames(express) == "ENSG00000065328", ]
ENSG00000065328values <- log2(ENSG00000065328values + 1)
ENSG00000065328 <- ifelse(ENSG00000065328values > median(ENSG00000065328values), 1, 0)
table(ENSG00000065328)
hr_samples <- sum(ENSG00000065328 == 1)
lr_samples <- sum(ENSG00000065328 == 0)

### Selecting Kidney survival object
KidneySurv <- survivalDOD

### Fitting survival model according ENSG00000065328
kidneyRNASurvENSG00000065328 <- survfit(KidneySurv ~ ENSG00000065328)

### Drawing Kaplan-Meier curve
jpeg(filename = "../figs/Kidney_KM_curve_ENSG00000065328_italicized_offname.jpg", width = 4000, height = 4000, res = 300)
ggsurvkm_ENSG00000065328 <- ggsurvplot(kidneyRNASurvENSG00000065328, data = KidneySurv, 
                       pval = TRUE, pval.method = TRUE,
                       xlab = "", ylab = "Survival probability", pval.size = 10,
                       #title = "Gene ENSG00000065328",
                       title = "MCM10",
                       legend.lab = c(paste0("Low-expression, n = ", hr_samples), 
                                      paste0("High-expression, n = ", lr_samples)), 
                       legend = c(0.75,0.8),
                       palette = c("skyblue1", "deeppink1"), size = 2,
                       #font.x = 5, font.y = 5, 
                       pval.coord = c(1000 ,0.1), 
                       pval.method.coord = c(0, 0.1),
                       risk.table = TRUE,
                       risk.table.height = 0.2,
                       tables.y.text = FALSE,
                       fontsize = 10,
                       )
ggsurvkm_ENSG00000065328 <- ggpar(ggsurvkm_ENSG00000065328,
                                  font.title    = c(34, "plain", "black"),         
                                  font.caption  = c(14, "plain"),        
                                  font.x        = c(30, "italic"),          
                                  font.y        = c(30, "italic"),      
                                  font.xtickslab = c(30, "plain"),
                                  font.ytickslab = c(30, "plain"),
                                  font.legend = c(30, "plain", "black"),
                                  #legend = "top"
)
ggsurvkm_ENSG00000065328$plot <- ggsurvkm_ENSG00000065328$plot + 
   #ggtitle("Gene ENSG00000065328") +
   ggtitle("MCM10") +
   theme(
     plot.title = element_text(face = "italic", color = "black", size=34, hjust = 0.5),
#     axis.text.x=element_text(color = "black", size=30),
#     axis.text.y=element_text(color = "black", size=30),
#     axis.title.x=element_text(size = 40),
#     axis.title.y=element_text(size = 40),
#     legend.text=element_text(size = 34),
#     legend.title=element_text(size = 34),
#     legend.key.size = unit(1.5, "cm")
)
ggsurvkm_ENSG00000065328
dev.off()

###########################################################################
### Creating prognostic variable for CATG00000107122

### loading elncRNA eSet and extracting expression matrix
load("~/DropboxMech/FANTOM6/TCGA/DGE/Kidney/objs/KidneyelncRNA.rda")
express <- exprs(KidneyelncRNA)
express <- express[, colnames(express) %in% sampleIDs]
### Consistency check: Expression colnames are in the same order as pheno rownames
all(rownames(pheno) == colnames(express))
rm(KidneyelncRNA)

### selecting expression
CATG00000107122values <- express[rownames(express) == "CATG00000107122", ]
CATG00000107122values <- log2(CATG00000107122values + 1)
CATG00000107122 <- ifelse(CATG00000107122values > median(CATG00000107122values), 1, 0)
table(CATG00000107122)
hr_samples <- sum(CATG00000107122 == 1)
lr_samples <- sum(CATG00000107122 == 0)

### Fitting survival model according CATG00000107122
kidneyRNASurvCATG00000107122 <- survfit(KidneySurv ~ CATG00000107122)

### Drawing Kaplan-Meier curve
jpeg(filename = "../figs/Kidney_KM_curve_CATG00000107122_italicized_offname.jpg", width = 4000, height = 4000, res = 300)
ggsurvkm_CATG00000107122 <- ggsurvplot(kidneyRNASurvCATG00000107122, data = KidneySurv, 
                                       pval = TRUE, pval.method = TRUE,
                                       xlab = "", ylab = "", pval.size = 10,
                                       title = "CATG00000107122",
                                       legend.lab = c(paste0("Low-expression, n = ", hr_samples), 
                                                      paste0("High-expression, n = ", lr_samples)), 
                                       legend = c(0.75,0.8),
                                       palette = c("skyblue1", "deeppink1"), size = 2,
                                       #font.x = 5, font.y = 5, 
                                       pval.coord = c(1000 ,0.1), 
                                       pval.method.coord = c(0, 0.1),
                                       risk.table = TRUE,
                                       risk.table.height = 0.2,
                                       tables.y.text = FALSE,
                                       fontsize = 10,
                                       )
ggsurvkm_CATG00000107122 <- ggpar(ggsurvkm_CATG00000107122,
                                  font.title    = c(34, "plain", "black"),         
                                  font.caption  = c(14, "plain"),        
                                  font.x        = c(30, "italic"),          
                                  font.y        = c(30, "italic"),      
                                  font.xtickslab = c(30, "plain"),
                                  font.ytickslab = c(30, "plain"),
                                  font.legend   = c(30, "plain", "black"),
                                  )
ggsurvkm_CATG00000107122$plot <- ggsurvkm_CATG00000107122$plot + 
  ggtitle("CATG00000107122") +
  theme(
    plot.title = element_text(face = "italic", color = "black", size=34, hjust = 0.5),
    #     axis.text.x=element_text(color = "black", size=30),
    #     axis.text.y=element_text(color = "black", size=30),
    #     axis.title.x=element_text(size = 40),
    #     axis.title.y=element_text(size = 40),
    #     legend.text=element_text(size = 34),
    #     legend.title=element_text(size = 34),
    #     legend.key.size = unit(1.5, "cm")
  )
ggsurvkm_CATG00000107122
dev.off()

###########################################################################
### Creating prognostic variable for ENSG00000235989

### loading elncRNA eSet and extracting expression matrix
load("~/DropboxMech/FANTOM6/TCGA/DGE/Kidney/objs/KidneydivergentP.rda")
express <- exprs(KidneydivergentP)
express <- express[, colnames(express) %in% sampleIDs]
### Consistency check: Expression colnames are in the same order as pheno rownames
all(rownames(pheno) == colnames(express))
rm(KidneydivergentP)

### selecting expression
ENSG00000235989values <- express[rownames(express) == "ENSG00000235989", ]
ENSG00000235989values <- log2(ENSG00000235989values + 1)
ENSG00000235989 <- ifelse(ENSG00000235989values > median(ENSG00000235989values), 1, 0)
table(ENSG00000235989)
hr_samples <- sum(ENSG00000235989 == 1)
lr_samples <- sum(ENSG00000235989 == 0)

### Fitting survival model according ENSG00000235989
kidneyRNASurvENSG00000235989 <- survfit(KidneySurv ~ ENSG00000235989)

### Drawing Kaplan-Meier curve
jpeg(filename = "../figs/Kidney_KM_curve_ENSG00000235989_italicized_offname.jpg", width = 4000, height = 4000, res = 300)
ggsurvkm_ENSG00000235989 <- ggsurvplot(kidneyRNASurvENSG00000235989, data = KidneySurv, 
                                       pval = TRUE, pval.method = TRUE,
                                       xlab = "Overall survival in days", ylab = "Survival probability", 
                                       pval.size = 10,
                                       #title = "Gene ENSG00000235989",
                                       title = "MORC2-AS",
                                       legend.lab = c(paste0("Low-expression, n = ", hr_samples), 
                                                      paste0("High-expression, n = ", lr_samples)), 
                                       legend = c(0.75,0.8),
                                       palette = c("skyblue1", "deeppink1"), size = 2,
                                       #font.x = 5, font.y = 5, 
                                       pval.coord = c(1000 ,0.1), 
                                       pval.method.coord = c(0, 0.1),
                                       risk.table = TRUE,
                                       risk.table.height = 0.2,
                                       tables.y.text = FALSE,
                                       fontsize = 10,
)
ggsurvkm_ENSG00000235989 <- ggpar(ggsurvkm_ENSG00000235989,
                                  font.title    = c(34, "plain", "black"),         
                                  font.caption  = c(14, "plain"),        
                                  font.x        = c(30, "italic"),          
                                  font.y        = c(30, "italic"),      
                                  font.xtickslab = c(30, "plain"),
                                  font.ytickslab = c(30, "plain"),
                                  font.legend   = c(30, "plain", "black"),
)
ggsurvkm_ENSG00000235989$plot <- ggsurvkm_ENSG00000235989$plot + 
  #ggtitle("Gene ENSG00000235989") +
  ggtitle("MORC2-AS") +
  theme(
    plot.title = element_text(face = "italic", color = "black", size=34, hjust = 0.5),
    #     axis.text.x=element_text(color = "black", size=30),
    #     axis.text.y=element_text(color = "black", size=30),
          axis.title.x=element_blank(),
    #     axis.title.y=element_text(size = 40),
    #     legend.text=element_text(size = 34),
    #     legend.title=element_text(size = 34),
    #     legend.key.size = unit(1.5, "cm")
  )
ggsurvkm_ENSG00000235989
dev.off()

###########################################################################
### Creating prognostic variable for ENSG00000196756

### loading elncRNA eSet and extracting expression matrix
load("~/DropboxMech/FANTOM6/TCGA/DGE/Kidney/objs/KidneyintergenicP.rda")
express <- exprs(KidneyintergenicP)
express <- express[, colnames(express) %in% sampleIDs]
### Consistency check: Expression colnames are in the same order as pheno rownames
all(rownames(pheno) == colnames(express))
rm(KidneyintergenicP)

### selecting expression
ENSG00000196756values <- express[rownames(express) == "ENSG00000196756", ]
ENSG00000196756values <- log2(ENSG00000196756values + 1)
ENSG00000196756 <- ifelse(ENSG00000196756values > median(ENSG00000196756values), 1, 0)
table(ENSG00000196756)
hr_samples <- sum(ENSG00000196756 == 1)
lr_samples <- sum(ENSG00000196756 == 0)

### Fitting survival model according ENSG00000196756
kidneyRNASurvENSG00000196756 <- survfit(KidneySurv ~ ENSG00000196756)

### Drawing Kaplan-Meier curve
jpeg(filename = "../figs/Kidney_KM_curve_ENSG00000196756_italicized_offname.jpg", width = 4000, height = 4000, res = 300)
ggsurvkm_ENSG00000196756 <- ggsurvplot(kidneyRNASurvENSG00000196756, data = KidneySurv, 
                                       pval = TRUE, pval.method = TRUE,
                                       xlab = "Overall survival in days", ylab = "", 
                                       pval.size = 10,
                                       #title = "Gene ENSG00000196756",
                                       title = "SNHG1",
                                       legend.lab = c(paste0("Low-expression, n = ", hr_samples), 
                                                      paste0("High-expression, n = ", lr_samples)), 
                                       legend = c(0.75,0.8),
                                       palette = c("skyblue1", "deeppink1"), size = 2,
                                       #font.x = 5, font.y = 5, 
                                       pval.coord = c(1000 ,0.1), 
                                       pval.method.coord = c(0, 0.1),
                                       risk.table = TRUE,
                                       risk.table.height = 0.2,
                                       tables.y.text = FALSE,
                                       fontsize = 10,
)
ggsurvkm_ENSG00000196756 <- ggpar(ggsurvkm_ENSG00000196756,
                                  font.title    = c(34, "plain", "black"),         
                                  font.caption  = c(14, "plain"),        
                                  font.x        = c(30, "italic"),          
                                  font.y        = c(30, "italic"),      
                                  font.xtickslab = c(30, "plain"),
                                  font.ytickslab = c(30, "plain"),
                                  font.legend   = c(30, "plain", "black"),
)
ggsurvkm_ENSG00000196756$plot <- ggsurvkm_ENSG00000196756$plot + 
  #ggtitle("Gene ENSG00000196756") +
  ggtitle("SNHG1") +
  theme(
    plot.title = element_text(face = "italic", color = "black", size=34, hjust = 0.5),
    #     axis.text.x=element_text(color = "black", size=30),
    #     axis.text.y=element_text(color = "black", size=30),
          axis.title.x=element_blank(),
    #     axis.title.y=element_text(size = 40),
    #     legend.text=element_text(size = 34),
    #     legend.title=element_text(size = 34),
    #     legend.key.size = unit(1.5, "cm")
  )
ggsurvkm_ENSG00000196756
dev.off()

