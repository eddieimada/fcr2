### Load required libraries
library(GenomicRanges)
library(SummarizedExperiment)
library(recount)

### Download data
download_study('GTEX', type = 'rse-fc')

### Scale data
GTEX <- scale_counts(rse_fc, by="auc")

### Extract counts matrix
rseCounts <- assays(GTEX)$counts
### Extract feature information
grList <- rowRanges(GTEX)
fdata <- mcols(grList)
### Creating list of pheno data splited by cancer type
pheno <- as.data.frame(colData(GTEX))

### Check number of samples and projects available
table(pheno$gdc_cases.project.project_id)

### Create expression set
gtexEset <- ExpressionSet(assayData = rseCounts,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(fdata))

save(gtexEset, file="gtexEset.rda")
