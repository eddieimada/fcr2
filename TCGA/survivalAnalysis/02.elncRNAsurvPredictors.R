###########################################################################
### Survival analysis by TCGA cancer types using e-lncRNA
### 
### Diego F Sanchez and Luigi Marchionni
###
###########################################################################
### First steps
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

###########################################################################
### Useful functions
###########################################################################

### Function to save data.frames from a list to files named as in the list names()
saveDFtoFiles <- function(x, y) {
  # x = names(list), y = list
  temp.df <- y[[x]]
  temp.df <- temp.df[ order(abs(temp.df$HR), decreasing=TRUE), ]
  temp.df <- temp.df[ order(temp.df$BH), ]
  write.csv(temp.df, file = paste("../text/SigEnhancerlncRNA_", x, ".csv", sep = ""))
}

###########################################################################
### Analysis by TCGA cancer type
###########################################################################

### Loading object with Cox HR
load("../objs/byOrganTCGATypesCoxelncRNA.rda")
load("../objs/byOrganTCGATypesSurvModels.rda")

###########################################################################
### Structure of the loaded objects:
class(coxModelByTypes)
names(coxModelByTypes)
str(coxModelByTypes[1])
head(coxModelByTypes[[1]])

class(survModelByTypes)
names(survModelByTypes)
str(survModelByTypes[1])
survModelByTypes[1]

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
# Saving into significantEnhancersByTCGAtypes.csv
write.csv(summaryDF, file = "../text/significantEnhancersByTCGAtypes.csv")


### Making plot for number of genes vs TCGA cancer types
jpeg(filename = "../figs/NumberOfEnhancersByTCGAtypes.jpg", width = 854)
ggplot(data = summaryDF, aes(x = Types, y = `FDR <= 0.05`)) +
  theme_classic()+
  geom_bar(stat = "identity", fill = "salmon") + xlab("") + ylab("") +
  theme(axis.text.x=element_text(color = "black", size=10, angle=90, vjust=.8, hjust=0.8), 
        plot.title = element_text(hjust=0.5))+
  ggtitle("Number of significant elncRNAs by TCGA cancer types")
dev.off()

###########################################################################
### data.frames with Significant genes (BH <= 0.05)
significantEnhancers <- lapply(coxModelByTypes, function(x){
  tmp <- as.data.frame(x[!is.na(x[, "HR"]), ]) # It is generating NAs because there are NAs in the Cox Models.
  tmp <- tmp[tmp[,"BH"] <= 0.05, ]
  })

### Saving data.frames with elncRNAs by TCGA types
invisible(lapply(names(significantEnhancers), saveDFtoFiles, y = significantEnhancers))

###########################################################################
### Common genes among subtypes
# Retaining just names
elncNames <- lapply(significantEnhancers, row.names)
elncNames <- elncNames[sapply(elncNames, length) > 0]

# elncRNAs shared among all cancer types
elncNamesCommon <- Reduce(f="intersect", x = elncNames)
elncNamesCommon

# Predictive elncRNAs
elncNamesPredictives <- Reduce(f="union", x = elncNames)
length(elncNamesPredictives)
write.csv(x = elncNamesPredictives, file = "../text/union_of_predictive_enhancers_through_tumors_list.csv")

# Comparing Predictive elncRNAs with DGE lncRNAs
DGElncRNAs <- c("CATG00000107122", "CATG00000054627", "CATG00000039286", 
                "ENSG00000231246", "ENSG00000255958")
DGElncRNAs %in% elncNamesPredictives

# Conparing Predictive with DGE lncRNAs by cancer type
lapply(significantEnhancers, function(x){DGElncRNAs %in% row.names(x)})

pairs <- combinations(length(names(elncNames)), 2, v = names(elncNames))

sharedEnhancersByPairs <- list()
for (i in 1:nrow(pairs)){
  sharedEnhancersByPairs[[ paste(pairs[i,1], pairs[i,2], sep = "-") ]] <- intersect(
    elncNames[[which(names(elncNames) == pairs[i,1])]], 
    elncNames[[which(names(elncNames) == pairs[i,2])]])
  print(paste(i, pairs[i,1], pairs[i,2], sep = " "))
}

countEnhancersByTypePairs <- as.data.frame(sapply(sharedEnhancersByPairs, length))
colnames(countEnhancersByTypePairs) <- "Predictive elncRNAs shared"
countEnhancersByTypePairs <- cbind(`Enhancer pairs` = row.names(countEnhancersByTypePairs), countEnhancersByTypePairs)
countEnhancersByTypePairs <- countEnhancersByTypePairs[order(countEnhancersByTypePairs$`Predictive elncRNAs`, 
                                                             decreasing = TRUE), ]
`Enhancer pairs` = factor(row.names(countEnhancersByTypePairs), levels = row.names(countEnhancersByTypePairs))
countEnhancersByTypePairs$`Enhancer pairs` <- `Enhancer pairs`
row.names(countEnhancersByTypePairs)<- 1:nrow(countEnhancersByTypePairs)

### Saving data.frame with elncRNAs by paired TCGA types
write.csv(countEnhancersByTypePairs, file = "../text/elncRNApairedBySubtypes.csv")

###########################################################################
### Analysis corresponding to enhancer 22: ENSG00000272666
###########################################################################

### identifying presence and significance level of enhancer 22.
sapply(coxModelByTypes, function(x){
  x[which("ENSG00000272666" == row.names(x)), ]
  })
sapply(significantEnhancers, function(x){"ENSG00000272666" %in% row.names(x)})

### ONLY KIDNEY DATA LOOK PROMISING. 
### Preparing Kidney cancer data
### Loading Kidney data
load("/Users/diegosanchez/DropboxMech/FANTOM6/TCGA/DGE//Kidney/objs/KidneyelncRNA.rda")
### extracting pheno data
pheno <- pData(KidneyelncRNA)
### keeping complete columns in pheno data
keep.cols <- which(apply(pheno, 2, function(x) !all(is.na(x))))
pheno <- pheno[,keep.cols]
### extracting expression data
express <- exprs(KidneyelncRNA)
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

### Creating prognostic variable for ENSG00000272666
ENSG00000272666values <- express[rownames(express) == "ENSG00000272666", ]
ENSG00000272666 <- ifelse(ENSG00000272666values > median(ENSG00000272666values), 1, 0)

### Selecting Kidney survival object
load("../objs/byOrganTCGATypesSurvObjects.rda")
KidneySurv <- survObjectByTypes$Kidney
rm(survObjectByTypes)

### Fitting survival model according ENSG00000272666
kidneyEnhancerSurv <- survfit(KidneySurv ~ ENSG00000272666)

### Drawing Kaplan-Meier curve
jpeg(filename = "../figs/KidneyKMcurve.jpg", width = 4000, height = 4000, res = 300)
ggsurvkm <- ggsurvplot(kidneyEnhancerSurv, data = KidneySurv, pval = TRUE, pval.method = TRUE,
           xlab = "Overall survival in days", ylab = "Survival probability", pval.size = 10,
           legend.lab = c("Low-risk, n = 442", "High-risk, n = 439"), legend = c(0.8,0.8),
           palette = c("blue", "red"), size = 2,
           #font.x = 5, font.y = 5, 
           pval.coord = c(1500 ,0.1), pval.method.coord = c(500, 0.1))
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

###########################################################################
### closing and saving session info
###########################################################################

### Session date and info
date()
sessionInfo()

### Clean
rm(list=ls())
quit(save = "no")
           


