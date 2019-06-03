#setwd("/users/eimada/FANTOM6/")

### Load libraries
library(SummarizedExperiment)
library(recount)

### Load objs
load("objs/rse_TCGA_pos.Rdata")
load("objs/uniqRangesPos.rda")
load("objs/fDataCAT.rda")

### Scale library to 40M
rse_scaled <- scale_counts(rse, round = FALSE)
tcga <- rse_scaled
rm(rse)
### Extract counts matrix
rseCounts <- assays(rse_scaled)$counts
rownames(rseCounts) <- rowRanges(rse_scaled)$ID

### Subset unambiguous ranges
rseCounts <- rseCounts[as.vector(uniqRangesPos$ranges),]

### Check order
all(rownames(rseCounts) == as.vector(uniqRangesPos$ranges))

### Summarize counts at gene level
nms <- colnames(rseCounts)
rseCountsGene <- apply(rseCounts, 2, function(x){
        tapply(x, uniqRangesPos$genes, FUN=sum)
        })
names(rseCountsGene) <- nms



### Get phenotypes
pheno <- as.data.frame(colData(tcga))
#pheno <- as.data.frame(lapply(pheno, unlist), stringsAsFactors=FALSE)

### Add rownames
#rownames(pheno) <- pheno$run

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
smp <- intersect(colnames(rseCountsGene), rownames(pheno))
str(smp)

### Available genes (common between expression and annotation table)
ids <- intersect(rownames(rseCountsGene), ann$geneID)
str(ids)

### Reorder and subset
pheno <- pheno[smp, ]
ann <- ann[ ids, ]
rseCountsGene <- rseCountsGene[ ids , smp ]

### Check order of samples
all(colnames(rseCountsGene) %in% rownames(pheno))
all(colnames(rseCountsGene) == rownames(pheno))

### Check order of genes
all(rownames(rseCountsGene) %in% rownames(ann))
all(rownames(rseCountsGene) == rownames(ann))


#################################################################
### Prepare MIAME
miame <- new("MIAME", name="GTEX_CAT",
             title="TCGA expression data - CAT FANTOM genes - Recount2",
             lab="Marchionni-Langmead-Leek-Jaffe")


#################################################################
### Create the expressionSet
TCGAPos <- ExpressionSet(assayData <- rseCountsGene,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(ann),
                          MIAME = miame)
annotation(TCGAPos) <- "CAT.FANTOM"

### Check validity
validObject(TCGAPos)


#################################################################
### ### Save
### save(gtexEset, file="objs/gtexEset.rda")
save(TCGAPos, file="objs/TCGAPos.rda")


#################################################################
### Session information and clean quit

### Session
sessionInfo()

### Clean
rm(list=ls())
