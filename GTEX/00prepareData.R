### Load required libraries
library(GenomicRanges)
library(SummarizedExperiment)
library(recount)

### Download data
download_study('GTEX', type = 'rse-fc')

### Scale data
GTEX_raw <- rse_fc
GTEX <- scale_counts(rse_fc, by="auc")

### Extract counts matrix
rseCounts_raw <- assays(GTEX_raw)$counts
rseCounts <- assays(GTEX)$counts
### Extract feature information
grList <- rowRanges(GTEX)
fdata <- mcols(grList)
### Creating list of pheno data splited by cancer type
pheno <- as.data.frame(colData(GTEX))

### Check number of samples and projects available
table(pheno$gdc_cases.project.project_id)

### Create expressions sets
gtexEset <- ExpressionSet(assayData = rseCounts,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(fdata))
gtexEset_raw <- ExpressionSet(assayData = rseCounts_raw,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(fdata))

save(gtexEset, file="gtexEset.rda")
save(gtexEset, file="gtexEset_raw.rda")
