#################################################################
### Boxplot of the expression level and speciificty

#################################################################
### Clean and set wd
rm(list=ls())

setwd("/users/eimada/FANTOM6/")

#################################################################
### Load packages
require(Biobase)
require(entropy)

### Load Data 
# Path to GTEx object
path <- "~/FANTOM6/GTEX/gtexEset.rda"
load(path)
load(file = "~/FANTOM6/GTEX/gtexMaxAllLog.rda")
load(file = "~/FANTOM6/GTEX/gtexEntropyAll.rda")

#################################################################
### Create Tissue Facets 
TissueTypes <- factor(gtexEset$smtsd)
#################################################################
### Feature Types 
### Create lncRNA categories (as factors)
FeatureTypes <- factor(paste(featureData(gtexEset)$CAT_geneClass,featureData(gtexEset)$CAT_DHS_type))
levels(FeatureTypes)

### Create intergenic p-lncRNA
levels(FeatureTypes)[grep("lncRNA_intergenic DHS_promoter", levels(FeatureTypes))] <- "intergenic p-lncRNA"

###Create divergent p-lncRNA 
levels(FeatureTypes)[grep("lncRNA_divergent DHS_promoter", levels(FeatureTypes))] <- "divergent p-lncRNA"

### Create elncRNA
levels(FeatureTypes)[grep("lncRNA_divergent DHS_enhancer", levels(FeatureTypes))] <- "elncRNA" 
levels(FeatureTypes)[grep("lncRNA_intergenic DHS_enhancer", levels(FeatureTypes))] <- "elncRNA" 
levels(FeatureTypes)[grep("lncRNA_antisense DHS_enhancer", levels(FeatureTypes))] <- "elncRNA" 
levels(FeatureTypes)[grep("lncRNA_sense_intronic DHS_enhancer", levels(FeatureTypes))] <- "elncRNA" 

### Create other-lncrna category 
levels(FeatureTypes)[grep ("lncRNA_divergent DHS_dyadic", levels(FeatureTypes))] <- "other-lncRNA"
levels(FeatureTypes)[grep ("lncRNA_intergenic DHS_dyadic", levels(FeatureTypes))] <- "other-lncRNA"
levels(FeatureTypes)[grep ("lncRNA_antisense DHS_dyadic", levels(FeatureTypes))] <- "other-lncRNA"
levels(FeatureTypes)[grep ("lncRNA_sense_intronic DHS_dyadic", levels(FeatureTypes))] <- "other-lncRNA"

levels(FeatureTypes)[grep ("lncRNA_divergent not_DHS", levels(FeatureTypes))] <- "other-lncRNA" 
levels(FeatureTypes)[grep ("lncRNA_intergenic not_DHS", levels(FeatureTypes))] <- "other-lncRNA"
levels(FeatureTypes)[grep ("lncRNA_antisense not_DHS", levels(FeatureTypes))] <- "other-lncRNA"
levels(FeatureTypes)[grep ("lncRNA_sense_intronic not_DHS", levels(FeatureTypes))] <- "other-lncRNA"

levels(FeatureTypes)[grep ("lncRNA_antisense DHS_promoter", levels(FeatureTypes))] <- "other-lncRNA"
levels(FeatureTypes)[grep ("lncRNA_sense_intronic DHS_promoter", levels(FeatureTypes))] <- "other-lncRNA"

### Simplify other categories
levels(FeatureTypes)[grep ("coding_mRNA ", levels(FeatureTypes))] <- "codingmRNA"
levels(FeatureTypes)[grep ("NA ", levels(FeatureTypes))] <- "NA"
levels(FeatureTypes)[grep ("pseudogene ", levels(FeatureTypes))] <- "pseudogene"
levels(FeatureTypes)[grep ("sense ", levels(FeatureTypes))] <- "sense-overlap"
levels(FeatureTypes)[grep ("short_ncRNA ", levels(FeatureTypes))] <- "short-ncRNA"
levels(FeatureTypes)[grep ("small_RNA ", levels(FeatureTypes))] <- "smallRNA"
levels(FeatureTypes)[grep ("structural", levels(FeatureTypes))] <- "structuralRNA"
levels(FeatureTypes)[grep ("uncertain", levels(FeatureTypes))] <- "uncertain-coding"

#################################################################
### Save plots
png("./figs/ExpressionSpecificitylncRNAtypes.png",
    height = 1180, width = 5000, res = 300)

### Create panel to fit 4 plots 
layout(matrix(1:4, nrow=1))

#######################
### Make plot for codingmRNA and add genes to plot
sel <- FeatureTypes == "codingmRNA"

smoothScatter(gtexEntropyAll[sel], gtexMaxAllLog[sel],
              colramp=colorRampPalette(c("white", "red3")), 
              xlab = "expression specificity", nrpoints=0,
              ylab = "max expression across facets (log2, tpm)", main = "Coding mRNA", ylim=c(-5,20))

gns <- c("ACTB", "CD74", "IL7", "SOX1", "NANOG")
selGns <- which(featureData(gtexEset)@data$HGNC_symbol %in% gns)
gnsSymb <- featureData(gtexEset)@data$HGNC_symbol[selGns]

text(gtexEntropyAll[selGns], gtexMaxAllLog[selGns], 
     labels = gnsSymb)

### Add median lines for entropy
mCodingent <- median(gtexEntropyAll[sel], na.rm = T)
abline(v= mCodingent, col="red3", lty=2)

mCodingexp <- median(gtexMaxAllLog[sel], na.rm = T)
abline(h= mCodingexp, col="red3", lty=2)

### Add grid
grid(col = "gray60")

#################################################################
### Make plot for divergent p-lncRNA and add genes to plot
sel <- FeatureTypes == "divergent p-lncRNA"

smoothScatter(gtexEntropyAll[sel], gtexMaxAllLog[sel],
              colramp=colorRampPalette(c("white", "purple")), 
              xlab = "expression specificity", nrpoints=0,
              ylab = "max expression across facets (log2, tpm)", main = "Divergent p-lncRNA", ylim=c(-5,20))

gns <- c("ZFAS1", "TUG1", "FENDRR", "PCAT7", "UCHL1-AS1")
selGns <- which(featureData(gtexEset)@data$HGNC_symbol %in% gns)
gnsSymb <- featureData(gtexEset)@data$HGNC_symbol[selGns]

text(gtexEntropyAll[selGns], gtexMaxAllLog[selGns], 
     labels = gnsSymb)

### Add median Line
mDPent <- median(gtexEntropyAll[sel], na.rm = T)
abline(v= mDPent, col="purple", lty=2)

mDPexp <- median(gtexMaxAllLog[sel], na.rm = T)
abline(h= mDPexp, col="purple", lty=2)

### Add grid
grid(col = "gray60")

#################################################################
### Make plot for intergenic p-lncRNA and add genes to plot      
sel <- FeatureTypes == "intergenic p-lncRNA"
smoothScatter(gtexEntropyAll[sel], gtexMaxAllLog[sel],
              colramp=colorRampPalette(c("white", "cyan4")),
              xlab = "expression specificity",  nrpoints=0,
              ylab = "max expression across facets (log2, tpm)", main = "Intergenic p-lncRNA", ylim=c(-5,20))

gns <- c("MALAT1", "NEAT1", "XIST", "HAR1A", "BCAR4")
selGns <- which(featureData(gtexEset)@data$HGNC_symbol %in% gns)
gnsSymb <- featureData(gtexEset)@data$HGNC_symbol[selGns]

text(gtexEntropyAll[selGns], gtexMaxAllLog[selGns], 
     labels = gnsSymb)

### Add median line
mIPent <- median(gtexEntropyAll[sel], na.rm = T)
abline(v= mIPent, col="cyan4", lty=2)

mIPexp <- median(gtexMaxAllLog[sel], na.rm = T)
abline(h= mIPexp, col="cyan4", lty=2)

### Add grid
grid(col = "gray60")

#################################################################
### Make plot for elncRNA and add genes to plot
sel <- FeatureTypes == "elncRNA"

smoothScatter(gtexEntropyAll[sel], gtexMaxAllLog[sel],
              colramp=colorRampPalette(c("white", "green4")), 
              xlab = "expression specificity",  nrpoints=0,
              ylab = "max expression across facets (log2, tpm)", main = "elncRNA", ylim=c(-5,20))

gns <- c("EGOT", "LUCAT1", "HCCAT5", "TRERNA1", "DLX6-AS1")
selGns <- which(featureData(gtexEset)@data$HGNC_symbol %in% gns)
gnsSymb <- featureData(gtexEset)@data$HGNC_symbol[selGns]

text(gtexEntropyAll[selGns], gtexMaxAllLog[selGns], 
     labels = gnsSymb)

### Add median Line
mElent <- median(gtexEntropyAll[sel], na.rm = T)
abline(v= mElent, col="green4", lty=2)

mElexp <- median(gtexMaxAllLog[sel], na.rm = T)
abline(h= mElexp, col="green4", lty=2)

### Add grid
grid(col = "gray60")

### Dev off
dev.off()

#################################################################
### Create boxplot for lncRNA's

### Save plot
jpeg("./figs/boxplot of specificity and expression of lncRNA types.jpeg",
     height = 1000, width = 2000)

### Layout 
layout(matrix(c(1,2), ncol=2))

### Remove uncessary levels and reorder
FeatureTypes2 <- FeatureTypes
levels(FeatureTypes2)[c(2,6,7,8)] <- NA
levs <- c("codingmRNA", "divergent p-lncRNA", "intergenic p-lncRNA","elncRNA")
FeatureTypes2 <- factor(FeatureTypes2, levels = levs)

### Create boxplot for Entropy measure
boxplot(gtexEntropyAll~FeatureTypes2, col = c("red3", "purple", "cyan4", "green4"),
        main = "Specificity for 4 lncRNA groups", ylab = "expression specificity (log2, cpm")

### Add Grid
grid(col = "gray60")

### Create boxplot for Max expression measure 
gtexMaxAllLog[gtexMaxAllLog == -Inf] <-  NA
boxplot(gtexMaxAllLog~FeatureTypes2, col = c("red3", "purple", "cyan4", "green4"),
        main = "Maximum expression for 4 lncRNA groups", ylab = "max expression (log2, cpm)")

### Add grid
grid(col = "gray 60")

### Dev off
dev.off()

#################################################################
### Session information and clean, quit

### Session
sessionInfo()

### Clean
rm(list=ls())

### Quit
q("no")


