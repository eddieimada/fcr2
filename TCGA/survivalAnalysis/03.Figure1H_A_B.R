###########################################################################
### Relative Enhancer Expression Change
### Figure 1H of the Pan-Cancer Analysis of Enhancer
### Eddie L Imada & Diego F Sanchez
###
###########################################################################
### First steps
### Set working directory
setwd(".")

### Clean enviroments
rm(list=ls())

### Load libraries
library(Biobase)
require(ggplot2)


###########################################################################
### Useful functions


###########################################################################
### loadind esets and making enhancer expression dataframe by types
###########################################################################

files <- list.files("../../DGE/", recursive = TRUE,
                   pattern=".+elnc.+.rda$", full.names = TRUE)

### removing repeated objects
files <- files[-grep("DGElist", files)]
files <- files[-grep("Kidney/", files)]
files <- files[-grep("Lung/", files)]
files <- files[-grep("Uterus/", files)]
files <- files[-grep("Colon/", files)]
files <- files[-grep("Rectum", files)]
files <- files[-grep("UCS", files)]

### Expression data.frame for 18 cancer types
allEnhExpressiondf <- data.frame(tcgaID=character(),
                 subtype=character(), 
                 tumorRPM=numeric(), 
                 normalRPM=numeric())

### Expression data.frame for cancer types common with Chen et al.
commonEnhExpressiondf <- data.frame(tcgaIDCommon=character(),
                                    subtypeCommon=character(), 
                                    tumorRPMCommon=numeric(), 
                                    normalRPMCommon=numeric())

###########################################################################
### loadind esets and making enhancer expression dataframe by subtypes
###########################################################################

for(file in files) {
  ### indivivual files for test
  #file <- "/Users/diegosanchez/DropboxMech/FANTOM6/TCGA/DGE//Stomach/objs/StomachelncRNA.rda"
  #file <- "/Users/diegosanchez/DropboxMech/FANTOM6/TCGA/DGE//Bladder/objs/BladderelncRNA.rda" #duplicated primary
  #file <- "/Users/diegosanchez/DropboxMech/FANTOM6/TCGA/DGE//Uterus_sep/objs/UCSelncRNA.rda" #not normal samples
  #file <- "/Users/diegosanchez/DropboxMech/FANTOM6/TCGA/DGE//Colorectal/objs/ColorectalelncRNA.rda"
  
  load(file)
  type <- gsub(".+/", "", file)
  type <- gsub(".rda", "", type)
  eset <- get(type)
  
  ###########################################################################
  ### extracting pheno and expression data
  pheno <- pData(eset)
  express <- exprs(eset)
  
  ### keeping complete columns in pheno data
  keep.cols <- which(apply(pheno, 2, function(x) !all(is.na(x)) ))
  pheno <- pheno[,keep.cols]
  
  ### available sample types (normal, primary, mts, etc)
  table(pheno$gdc_cases.samples.sample_type)
  
  ##########################################################################
  ### keeping normal and primary tumor samples in Pheno data
  pheno <- pheno[which(pheno$gdc_cases.samples.sample_type == "Solid Tissue Normal" | 
                   pheno$gdc_cases.samples.sample_type == "Primary Tumor"), ]
  
  ## Extracting normal samples TCGA ID
  tcgaSampleID <- pheno[which(pheno$gdc_cases.samples.sample_type == 
                                "Solid Tissue Normal"),]$gdc_cases.submitter_id
  ## Subsetting pheno data with same TCGA ID as normal samples
  keepTcgaSampleID <- which(pheno$gdc_cases.submitter_id %in% tcgaSampleID)
  pheno <- pheno[keepTcgaSampleID,]
  tmpTable <- table(pheno$gdc_cases.submitter_id, pheno$gdc_cases.samples.sample_type)
  pheno <- pheno[which(pheno$gdc_cases.submitter_id 
              %in% rownames(tmpTable[(tmpTable[,1] & tmpTable[,2]), ])), ]
  
  ## available sample types based on normal TCGA IDs
  table(pheno$gdc_cases.samples.sample_type)
  ## sample with frozen and ffpe could be available, delete ffpe
  idx <- which(pheno$gdc_cases.samples.is_ffpe == TRUE)
  if (length(idx) != 0){
    pheno <- pheno[-idx,]
  }
  
  ## looking for cases with duplicated id (probable normal vs primary tumor)
  dup <- unique(pheno$gdc_cases.submitter_id[duplicated(pheno$gdc_cases.submitter_id)])
  
  ## Some cases have more than 1 primary tumor
  dupTable <- table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup])
  remove <- names(dupTable)[dupTable > 2]
  
  ## this case DELETE ALL duplicates (NOT implemented)
  # dup <- dup[!dup %in% remove]
  
  ## this case pick one occurence chosen by random
  if(length(remove) > 0){
  for (removing in remove) {
    x <- which(pheno$gdc_cases.submitter_id %in% removing & 
                 pheno$gdc_cases.samples.sample_type == "Primary Tumor")
    set.seed(1234)
    toDrop <- sample(x, length(x) -1)
    pheno <- pheno[-toDrop, ]
  }
  }
  table(pheno$gdc_cases.submitter_id[pheno$gdc_cases.submitter_id %in% dup] ) 
  
  ## available sample types based on normal TCGA IDs (Matching normal vs tumor)
  table(pheno$gdc_cases.samples.sample_type) 
  # check if quantity of normal == tumor
  normalToPrimaryTab <- table(pheno$gdc_cases.submitter_id, pheno$gdc_cases.samples.sample_type) 
  normalToPrimaryTab
  # ordering by sample ID and sample type (tumor vs normal) 
  pheno <- pheno[with(pheno, order(gdc_cases.submitter_id, gdc_cases.samples.sample_type)), ]
  
  ### Keeping normal and primary tumor samples in Expression data
  sampleIDs <- rownames(pheno)
  express <- express[, colnames(express) %in% sampleIDs]
  ### Expression colnames are in the same order as pheno rownames
  all(rownames(pheno) == colnames(express))
  express<- express[,row.names(pheno)]
  # recheck
  all(rownames(pheno) == colnames(express))
  
  ### CREATING THE GLOBAL EXPRESSION LEVEL FOR ALL THE ENHANCERS
  ### Mean of RPM by sample
  
  ## set for deletion of uninformative enhancers ## not done
  # uninformative <- which(apply(express, 1, function(x){sum(x) == 0}))
  # length(uninformative)
  # deleting uninformative enhancers
  # express <- express[-uninformative, ]
  
  # obtaining mean
  expressMean <- apply(express, 2, mean)
  
  ### Creating and binding data.frame of GLOBAL EXPRESSION LEVELS by all types
  indexNorm <- which(pheno$gdc_cases.samples.sample_type == "Solid Tissue Normal")
  indexTumor <- which(pheno$gdc_cases.samples.sample_type == "Primary Tumor")
  # checking coorespondence between normal and tu
  all(rownames(pheno[indexNorm,]) == names(expressMean[indexNorm]))
  all(rownames(pheno[indexTumor,]) == names(expressMean[indexTumor]))
  tumorRPM <- expressMean[indexTumor]
  normalRPM <- expressMean[indexNorm]
  tcgaID <- unique(pheno$gdc_cases.submitter_id)
  subtype <- gsub(".+-", "", pheno$gdc_cases.project.project_id)
  subtype <- subtype[1:length(tcgaID)]
  df <- cbind.data.frame(tcgaID, subtype, tumorRPM, normalRPM)
  allEnhExpressiondf <- rbind(df, allEnhExpressiondf)
  
  ### CREATING THE GLOBAL EXPRESSION LEVEL FOR COMMON SET OF ENHANCERS
  ### Subsetting by common set.
  mappingTable <- read.csv("../../../Manuscript/Pan_cancer/text/mapping_table_andersson.csv", header = T,
                           stringsAsFactors = FALSE)
  commonEnh <- unique(mappingTable$fantom_cat)
  expressCommon <- express[rownames(express) %in% as.character(commonEnh),]
  
  ### Mean of RPM by sample
  expressMeanCommon <- apply(expressCommon, 2, mean)
  
  ### Creating and binding the data.frame
  indexNormCommon <- which(pheno$gdc_cases.samples.sample_type == "Solid Tissue Normal")
  indexTumorCommon <- which(pheno$gdc_cases.samples.sample_type == "Primary Tumor") 
  all(rownames(pheno[indexNorm,]) == names(expressMeanCommon[indexNorm]))
  all(rownames(pheno[indexTumor,]) == names(expressMeanCommon[indexTumor]))
  tumorRPMCommon <- expressMeanCommon[indexTumor]
  normalRPMCommon <- expressMeanCommon[indexNorm]
  tcgaIDCommon <- unique(pheno$gdc_cases.submitter_id)
  subtypeCommon <- gsub(".+-", "", pheno$gdc_cases.project.project_id)
  subtypeCommon <- subtypeCommon[1:length(tcgaIDCommon)]
  df <- cbind.data.frame(tcgaIDCommon, subtypeCommon, tumorRPMCommon, normalRPMCommon)
  commonEnhExpressiondf <- rbind(df, commonEnhExpressiondf)
  
  ### Removing eset object
  rm(list = type) #use list to remove the value of type and no type
}

###########################################################################
### Final GLOBAL EXPRESSION data frame with all the subtypes
###########################################################################
write.csv(allEnhExpressiondf, file = "../text/globalAllEnhancerExpression.csv", sep = ",")
### Cancer subtypes
table(allEnhExpressiondf$subtype)
###########################################
### analysis by subtype in all enhancers
###########################################
subtypes <- names(table(allEnhExpressiondf$subtype))
### initializing data.frame
relativeChangeAllEnhdf <- data.frame(subtype = character(),
                        tumor = numeric(),
                        normal = numeric(),
                        relativeChange = numeric(),
                        p_value = numeric())
  
for (subtype in subtypes) {
    dftemp <- allEnhExpressiondf[which(allEnhExpressiondf$subtype == subtype), ]
    tumorMean <- mean(dftemp$tumorRPM)
    normalMean <- mean(dftemp$normalRPM)
    relativeChange <- ((tumorMean/normalMean) - 1)*100
    pairedtstat <- t.test(dftemp$tumorRPM, dftemp$normalRPM, paired=TRUE)
    tempList <- list("subtype" = subtype, "tumor" = tumorMean, "normal" = normalMean,
                 "relativeChange" = relativeChange, "p_value" = pairedtstat$p.value)
    relativeChangeAllEnhdf <- rbind(relativeChangeAllEnhdf, as.data.frame(tempList))
}

### GLOBAL ACTIVATION data.frame with all the types
write.csv(relativeChangeAllEnhdf, file = "../text/globalAllEnhancerActivation.csv", sep = ",")
  
###########################################
### analysis by subtype in common enhancers
###########################################
#subtypes <- names(table(allEnhExpressiondf$subtype))
relativeChangeCommonEnhdf <- data.frame(subtype = character(),
                                       tumor = numeric(),
                                       normal = numeric(),
                                       relativeChange = numeric(),
                                       p_value = numeric())
  
for (subtype in subtypes) {
  dftempCommon <- commonEnhExpressiondf[which(commonEnhExpressiondf$subtype == subtype), ]
  tumorMeanCommon <- mean(dftempCommon$tumorRPMCommon)
  normalMeanCommon <- mean(dftempCommon$normalRPMCommon)
  relativeChangeCommon <- ((tumorMeanCommon/normalMeanCommon) - 1)*100
  pairedtstatCommon <- t.test(dftempCommon$tumorRPMCommon, dftempCommon$normalRPMCommon, paired=TRUE)
  tempListCommon <- list("subtype" = subtype, "tumor" = tumorMeanCommon, "normal" = normalMeanCommon,
                     "relativeChange" = relativeChangeCommon, "p_value" = pairedtstatCommon$p.value)
  relativeChangeCommonEnhdf <- rbind(relativeChangeCommonEnhdf, as.data.frame(tempListCommon))
}
  
### SUBSETTING TO THE SAME CANCER TYPES AS IN CHEN et al.
drop <- c("CHOL", "UCEC", "KICH")
relativeChangeAllEnhdf <- relativeChangeAllEnhdf[-which(relativeChangeAllEnhdf$subtype %in% drop),]
relativeChangeCommonEnhdf <- relativeChangeCommonEnhdf[-which(relativeChangeCommonEnhdf$subtype %in% drop),]
  
### adjusting p_values
relativeChangeAllEnhdf <- cbind.data.frame(relativeChangeAllEnhdf, 
                                       adj_p_values = p.adjust(relativeChangeAllEnhdf$p_value, 
                                                               method = "BH"))
  
relativeChangeCommonEnhdf <- cbind.data.frame(relativeChangeCommonEnhdf, 
                                       adj_p_values = p.adjust(relativeChangeCommonEnhdf$p_value, 
                                                               method = "BH"))
  
### p_value significant to a certain level (0.05)
relativeChangeAllEnhdf <- cbind.data.frame(relativeChangeAllEnhdf, 
                                 p_sig = ifelse(relativeChangeAllEnhdf$adj_p_values > 0.05, 
                                                FALSE, TRUE))
relativeChangeCommonEnhdf <- cbind.data.frame(relativeChangeCommonEnhdf, 
                                p_sig = ifelse(relativeChangeCommonEnhdf$adj_p_values > 0.05, 
                                               FALSE, TRUE))  
  
### relative change is positive or negative (up- or down-regulated respectively)
relativeChangeAllEnhdf <- cbind.data.frame(relativeChangeAllEnhdf, 
                                up_regulated = ifelse(relativeChangeAllEnhdf$relativeChange > 0, 
                                                             TRUE, FALSE))
relativeChangeCommonEnhdf <- cbind.data.frame(relativeChangeCommonEnhdf, 
                                up_regulated = ifelse(relativeChangeCommonEnhdf$relativeChange > 0, 
                                                             TRUE, FALSE)) 

### Final GLOBAL ACTIVATION data.frame with Chen et al types and BH FDR.
write.csv(relativeChangeAllEnhdf, file = "../text/globalAllEnhancersActivationFDR.csv", sep = ",")
write.csv(relativeChangeCommonEnhdf, file = "../text/globalCommonEnhancersActivationFDR.csv", sep = ",")


### levels for coloring the figure according regulation and significance
relativeChangeAllEnhdf <- cbind.data.frame(relativeChangeAllEnhdf, 
                              colorLevels = ifelse(relativeChangeAllEnhdf$up_regulated, 1, 2))
relativeChangeAllEnhdf$colorLevels[which(!relativeChangeAllEnhdf$p_sig)] <- 3
relativeChangeAllEnhdf$colorLevels <- factor(relativeChangeAllEnhdf$colorLevels)
levels(relativeChangeAllEnhdf$colorLevels) <- c("up-regulated", "down-regulated", "FDR > 0.05")
relativeChangeAllEnhdf <- relativeChangeAllEnhdf[order(-relativeChangeAllEnhdf$relativeChange), ]
relativeChangeAllEnhdf$subtype <-factor(relativeChangeAllEnhdf$subtype, levels = relativeChangeAllEnhdf$subtype)

relativeChangeCommonEnhdf <- cbind.data.frame(relativeChangeCommonEnhdf, 
                              colorLevels = ifelse(relativeChangeCommonEnhdf$up_regulated, 1, 2))
relativeChangeCommonEnhdf$colorLevels[which(!relativeChangeCommonEnhdf$p_sig)] <- 3
relativeChangeCommonEnhdf$colorLevels <- factor(relativeChangeCommonEnhdf$colorLevels)
levels(relativeChangeCommonEnhdf$colorLevels) <- c("up-regulated", "down-regulated", "FDR > 0.05")
relativeChangeCommonEnhdf <- relativeChangeCommonEnhdf[order(-relativeChangeCommonEnhdf$relativeChange), ]
relativeChangeCommonEnhdf$subtype <-factor(relativeChangeCommonEnhdf$subtype, levels = relativeChangeCommonEnhdf$subtype)    
  
### Making plots
jpeg(filename = "../figs/allEnhancersActivationTCGAtypes.jpg", width = 854)
ggplot(data = relativeChangeAllEnhdf, aes(x = subtype, y = relativeChange, fill = colorLevels)) +
  geom_bar(stat="identity", width=0.7)+
  theme(axis.text.x=element_text(color = "black", size=10, angle=90, vjust=.8, hjust=0.8),
  legend.position= c(0.8,0.8), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Relative enhancer change using all the enhancers")+
  guides(fill=guide_legend(title=NULL, nrow = 3, byrow = TRUE))+
  labs(y = "Relative expression change", x = "TCGA cancer types")
dev.off()  

jpeg(filename = "../figs/commonEnhancersActivationTCGAtypes.jpg", width = 854)    
ggplot(data = relativeChangeCommonEnhdf, aes(x = subtype, y = relativeChange, fill = colorLevels)) +
  geom_bar(stat="identity", width=0.7)+
  theme(axis.text.x=element_text(color = "black", size=10, angle=90, vjust=.8, hjust=0.8),
        legend.position= c(0.8,0.8), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Relative enhancer change using common enhancers")+
  guides(fill=guide_legend(title=NULL, nrow = 3, byrow = TRUE))+
  labs(y = "Relative expression change", x = "TCGA cancer types")
dev.off()

###########################################################################
### Session date and info
date()
sessionInfo()

### Clean
rm(list=ls())
quit(save = "no")