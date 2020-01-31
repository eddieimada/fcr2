### Load required libraries
library(GenomicRanges)
library(SummarizedExperiment)
library(recount)

### Download data
download_study('TCGA', type = 'rse-fc')

### Scale data
TCGA <- scale_counts(rse_fc, by="auc")

### Extract counts matrix
rseCounts <- assays(TCGA)$counts
### Extract feature information
grList <- rowRanges(TCGA)
fdata <- mcols(grList)
### Creating list of pheno data splited by cancer type
pheno <- as.data.frame(colData(TCGA))

### Check number of samples and projects available
table(pheno$gdc_cases.project.project_id)

### Create expression set
TCGAeset <- ExpressionSet(assayData = rseCounts,
                      phenoData = AnnotatedDataFrame(pheno),
                      featureData= AnnotatedDataFrame(fdata))

save(TCGAeset, file="tcgaEset.rda")
