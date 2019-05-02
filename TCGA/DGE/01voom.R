###########################################################################
### Limma/Voom differential expression analysis correcting for batch effects
### Eddie Imada & Luigi Marchionni
### TCGA Prostate cancer divergent promoters
# This is an example script. This analysis applies to all cancers and RNA classes
###########################################################################
### Set working directory

### Clean enviroments
rm(list=ls())

### Load libraries
library(Biobase)
library(limma)
require(edgeR)
require(sva)
library(dplyr)
###########################################################################
### load eset
load("objs/ProstateCodingmRNA.rda")

###########################################################################
### Create design matrix
# Remove genes without information
keep <- complete.cases(fData(ProstateCodingmRNA)[,1:2])
fData(ProstateCodingmRNA) <- fData(ProstateCodingmRNA)[keep,]
# Remove metastatic samples
pheno <- pData(ProstateCodingmRNA)
table(pheno$gdc_cases.samples.sample_type)
metaIndex <- which(pheno$gdc_cases.samples.sample_type == "Metastatic")
pheno <- pheno[-metaIndex,]
express <- exprs(ProstateCodingmRNA)[,-metaIndex]
express <- express[keep,]

# Create factor and rename
groups <- pheno$gdc_cases.samples.sample_type
groups <- gsub("Solid Tissue Normal", "NORMAL", groups)
groups <- gsub("Primary Tumor", "TUMOR", groups)

# Relevel factor to be set as intercept
groups <- as.factor(groups)
groups <- relevel(groups,"NORMAL")

# Create null model for SVA analysis
dMat0 <- model.matrix( ~ 1 , data=pheno)

# Model matrix with intectept
dMat <- model.matrix( ~ groups)
colnames(dMat) <- c("Intercept", "TUMORvsNORMAL")


###########################################################################
###########################################################################
### GLM using voom() and limma
# Prepare expression matrix
dge <- as.matrix(round(express))

# Filter low counts genes
keepGns <- rowSums(dge > 5) >= 100
table(keepGns)
dge <- dge[keepGns, ]

# Create DGEList
dge <- DGEList(counts=dge, group=groups, genes=featureData(ProstateCodingmRNA)@data[keepGns,])
# TMM normalization by library size
dge <- calcNormFactors(dge, method="TMM")
# Note: observe voom trend, ideally it should decrease, if you see an upward trend
# at the beggining maybe do some more agressive filtering
de <- voom(dge, design = dMat, normalize.method = "none",
           plot = TRUE, span=0.4)


###########################################################################
###########################################################################
### SVA analysis
# Identify the number of surrogate variables
nSv <- num.sv(de$E, dMat, method="be")

# Estimate the surrogate variables
svobj <- sva(de$E, dMat, dMat0, n.sv=nSv)

# Final design (Here i'm only correcting for the first 3 SV, this is arbitrary)
# Adding SV to the design matrix
dMatSv <- cbind(dMat, svobj$sv[,1:3])
colnames(dMatSv)[colnames(dMatSv) %in% ""] <- paste("SV", 1:3, sep="")
colnames(dMatSv)[1] <- "Intercept"

###########################################################################
### With sva, using voomed data and limma: Fit generalized linear model
fit <- lmFit(de, dMatSv)
fit <- eBayes(fit)
plotSA(fit)

# Extract DEG genes table
mrnaPCA <- topTable(fit, coef= "TUMORvsNORMAL", n=Inf, genelist = fit$genes, resort.by="t")

###########################################################################
# Save DGE List and CSV
save(mrnaPCA, file="objs/DGE_mrna_prad.rda")
elncSignf <- filter(mrnaPCA, adj.P.Val < 0.01)
write.csv(elncSignf, file = paste0("./text/DGE_mrna_prad.csv"),
          row.names = FALSE)
###########################################################################
### Session information
sessionInfo()