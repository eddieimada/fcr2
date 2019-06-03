setwd("/users/eimada/FANTOM6/")

### Load libraries
library(SummarizedExperiment)
library(recount)

### Load objs
load("/dcl01/lmarchio/data/FANTOM6/CATdjn/djnUnstrand/DjnNeg/rse_SRP012682.Rdata")
load("objs/uniqRangesNeg.rda")
load("objs/fDataCAT.rda")

### Scale library to 40M
rse_scaled <- scale_counts(rse, round = FALSE)
gtex <- rse_scaled

### Extract counts matrix
rseCounts <- assays(rse_scaled)$counts
rownames(rseCounts) <- rowRanges(rse_scaled)$ID

### Subset unambiguous ranges
uniqRangesNeg <- uniqRangesNeg[as.vector(uniqRangesNeg$ranges) %in% rownames(rseCounts), ]
rseCounts <- rseCounts[as.vector(uniqRangesNeg$ranges),]

### Check order
all(rownames(rseCounts) == as.vector(uniqRangesNeg$ranges))

### Summarize counts at gene level
nms <- colnames(rseCounts)
rseCountsGene <- apply(rseCounts, 2, function(x){
        tapply(x, uniqRangesNeg$genes, FUN=sum)
        })
names(rseCountsGene) <- nms



### Get phenotypes
pheno <- as.data.frame(colData(gtex))
pheno <- as.data.frame(lapply(pheno, unlist), stringsAsFactors=FALSE)

### Add rownames
rownames(pheno) <- pheno$run

### Get rid of useless columns
keep <- sapply(pheno, function(x) length(unique(x)) > 1)
table(keep)

### Subset
pheno <- pheno[ , keep]


#################################################################
### Get gene annotation
ann <- fDataCAT[rownames(rseCountsGene),]

### Add rownames
rownames(ann) <- ann$geneID

### Get rid of useless columns
keep <- sapply(ann, function(x) length(unique(x)) > 1)
table(keep)

### Subset
ann <- ann[ , keep]


#################################################################
### ### Remove genes for which all samples have 0 counts
### rseCountsGene <- rseCountsGene[ rowSums(rseCountsGene) > 0 ,]

### Available samples (common between expression and phenotype table)
smp <- intersect(colnames(rseCountsGene), pheno$run)
str(smp)

### Available genes (common between expression and annotation table)
ids <- intersect(rownames(rseCountsGene), ann$geneID)
str(ids)

### Reorder and subset
pheno <- pheno[smp, ]
ann <- ann[ ids, ]
rseCountsGene <- rseCountsGene[ ids , smp ]

### Check order of samples
all(colnames(rseCountsGene) %in% pheno$run)
all(colnames(rseCountsGene) == pheno$run)

### Check order of genes
all(rownames(rseCountsGene) %in% rownames(ann))
all(rownames(rseCountsGene) == rownames(ann))


#################################################################
### Prepare MIAME
miame <- new("MIAME", name="GTEX_CAT",
             title="GTEX expression data - CAT FANTOM genes - Recount2",
             lab="Marchionni-Langmead-Leek-Jaffe")


#################################################################
### Create the expressionSet
gtexEsetNeg <- ExpressionSet(assayData <- rseCountsGene,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(ann),
                          MIAME = miame)
annotation(gtexEsetNeg) <- "CAT.FANTOM"

### Check validity
validObject(gtexEsetNeg)


#################################################################
### ### Save
### save(gtexEset, file="objs/gtexEset.rda")
save(gtexEsetNeg, file="objs/gtexEsetNeg.rda")


#################################################################
### Session information and clean quit

### Session
sessionInfo()

### Clean
rm(list=ls())
