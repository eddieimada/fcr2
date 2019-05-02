

#################################################################
###Tejasvi Matam
### Boxplot of the expression level and speciificty

#################################################################
### Clean and set wd
rm(list=ls())

#################################################################
### Load packages
require(Biobase)
require(entropy)

### Load Data 
# Path to GTEx object
path <- "~/FANTOM6/GTEX/gtexEset.rda"
load(path)
load("~/FANTOM6/GTEX/gtexTPM.rda")


### Load objs to run on laptop
#################################################################
### Create Tissue Facets 
TissueTypes <- factor(gtexEset$smtsd)

#################################################################
### Overall max
gtexMaxAll <- apply(gtexTPM,1, max)

### Log2
gtexMaxAllLog <- log2(gtexMaxAll)
save(gtexMaxAllLog, file = "objs/gtexMaxAllLogTPM.rda")
#################################################################
### Calculate Entropy
entropyFunc <- function(x, y) {
    sums <- tapply(x,y, mean, na.rm=TRUE)
    sums[sums < 0.1] <- 0
    sumsN <- entropy.empirical(sums, unit = "log2")/
        log2(length(levels(y)))
    out <- 1-sumsN
    out
}

# entropyFunc <- function(x, y) {
#     
#     sums <- tapply(x,y, mean, na.rm=TRUE)
#     sums <- ifelse(sums < 1, 0, 1)
#     #sums[sums < 1] <- 0
#     sumsN <- entropy.empirical(sums, unit = "log2")/
#         log2(length(levels(y)))
#     out <- 1-sumsN
# }



### Apply entropy function 
gtexEntropyTissue <- apply(gtexCPM, 1, entropyFunc, y= TissueTypes)
save(gtexEntropyTissue, file = "./objs/gtexEntropyTissues.rda")
### Apply empirical entropy
gtexEntropyAll <- -1*apply(gtexCPM, 1, entropy.empirical)

#################################################################
### Feature Types 
### Create lncRNA categories (as factors)
gtexfData <- fData(gtexEset)
FeatureTypes <- factor(paste(gtexfData$CAT_geneClass,gtexfData$CAT_DHS_type))
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
jpeg("./figs/ExpressionSpecificitylncRNAtypes.jpeg",
     height = 3000, width = 12000, res = 300)

### Create panel to fit 4 plots 
layout(mat = matrix(c(1,2,0,3,4,5,0,6,7,8,0,9,10,11,0,12), 2,8),  
       heights = c(2,8),
       widths = c(8,2.5,8,2.5,8,2.5,8,2.5))

#######################
### Make plot for codingmRNA and add genes to plot
sel <- FeatureTypes == "codingmRNA"
par(mar=c(0, 5, 1.1, 0),cex = 2.2)
corners <- par("usr")
boxplot(gtexEntropyTissue[sel], horizontal=TRUE , ylim=c(0,1), xaxt="n" , col
        ="firebrick1" , frame=F, outline = F, boxwex = 0.3, staplewex = 0.3,
        main = "Coding mRNA")
med <- median(gtexEntropyTissue[sel], na.rm = T)
text(x = med, y = corners[2]+0.3, labels = round(med,2), col = "red", srt = 0)
par(mar=c(5, 5, 0, 0))
smoothScatter(gtexEntropyTissue[sel], gtexMaxAllLog[sel],
              colramp=colorRampPalette(c("white", "firebrick1")),
              transformation = function(x) x,
              nrpoints = 0,
              ylim = c(-5,20), xlim = c(0,1),
              xlab = "expression specificity",
              ylab = "max expression across facets ")
grid (NULL,NULL, lty = 2, col = "grey60") 

gns <- c("ACTB", "CD74", "IL7", "SOX1", "NANOG")
selGns <- which(gtexfData$HGNC_symbol %in% gns)
gnsSymb <- gtexfData$HGNC_symbol[selGns]
points(gtexEntropyTissue[selGns], gtexMaxAllLog[selGns], 
       pch = 20, col = "black")
text(gtexEntropyTissue[selGns], gtexMaxAllLog[selGns], 
     labels = gnsSymb, pos = c(4,3,1,4), cex = 0.8)

### Add median lines for entropy
mCodingent <- median(gtexEntropyTissue[sel], na.rm = T)
abline(v= mCodingent, col="firebrick1", lty=2)
mCodingexp <- median(gtexMaxAllLog[sel], na.rm = T)
abline(h= mCodingexp, col="firebrick1", lty=2)

par(mar=c(5, 0, 0, 2))
corners <- par("usr")
boxplot(gtexMaxAllLog[sel], ylim=c(-5,20), yaxt="n" , col="firebrick1" , frame=F, outline = F, boxwex = 0.3, staplewex = 0.3)
med <- median(gtexMaxAllLog[sel], na.rm = T)
text(x = corners[2]+0.15, y = med, labels = round(med, 2), col = "firebrick1", srt = 270)




#################################################################
### Make plot for divergent p-lncRNA and add genes to plot
sel <- FeatureTypes == "divergent p-lncRNA"

par(mar=c(0, 5, 1.1, 0),cex = 2.2)
corners <- par("usr")
boxplot(gtexEntropyTissue[sel], horizontal=TRUE , ylim=c(0,1), xaxt="n" , col
        ="purple3" , frame=F, outline = F, boxwex = 0.3, staplewex = 0.3,
        main = "Divergent p-lncRNA")
med <- median(gtexEntropyTissue[sel], na.rm = T)
text(x = med, y = corners[2]-0.3, labels = round(med,2), col = "purple3", srt = 0)
par(mar=c(5, 5, 0, 0))
smoothScatter(gtexEntropyTissue[sel], gtexMaxAllLog[sel],
              colramp=colorRampPalette(c("white", "purple3")),
              transformation = function(x) x,
              nrpoints = 0,
              ylim = c(-5,20), xlim = c(0,1),
              xlab = "expression specificity",
              ylab = "max expression across facets ")
grid (NULL,NULL, lty = 2, col = "grey60") 

gns <- c("ZFAS1", "TUG1", "FENDRR", "PCAT7", "UCHL1-AS1")
selGns <- which(gtexfData$HGNC_symbol %in% gns)
gnsSymb <- gtexfData$HGNC_symbol[selGns]
points(gtexEntropyTissue[selGns], gtexMaxAllLog[selGns], 
       pch = 20, col = "black")
text(gtexEntropyTissue[selGns], gtexMaxAllLog[selGns], 
     labels = gnsSymb, pos = 4, cex = 0.8)

### Add median lines for entropy
mCodingent <- median(gtexEntropyTissue[sel], na.rm = T)
abline(v= mCodingent, col="purple3", lty=2)
mCodingexp <- median(gtexMaxAllLog[sel], na.rm = T)
abline(h= mCodingexp, col="purple3", lty=2)

par(mar=c(5, 0, 0, 2))
corners <- par("usr")
boxplot(gtexMaxAllLog[sel], ylim=c(-5,20), yaxt="n" , col="purple3" , frame=F, outline = F, boxwex = 0.3, staplewex = 0.3)
med <- median(gtexMaxAllLog[sel], na.rm = T)
text(x = corners[2]+0.15, y = med, labels = round(med, 2), col = "purple3", srt = 270)

#################################################################
### Make plot for intergenic p-lncRNA and add genes to plot      
sel <- FeatureTypes == "intergenic p-lncRNA"

par(mar=c(0, 5, 1.1, 0),cex = 2.2)
corners <- par("usr")
boxplot(gtexEntropyTissue[sel], horizontal=TRUE , ylim=c(0,1), xaxt="n" , col
        ="deepskyblue3" , frame=F, outline = F, boxwex = 0.3, staplewex = 0.3,
        main = "Intergenic p-lncRNA")
med <- median(gtexEntropyTissue[sel], na.rm = T)
text(x = med, y = corners[2]-0.3, labels = round(med,2), col = "deepskyblue3", srt = 0)
par(mar=c(5, 5, 0, 0))
smoothScatter(gtexEntropyTissue[sel], gtexMaxAllLog[sel],
              colramp=colorRampPalette(c("white", "deepskyblue3")),
              transformation = function(x) x,
              nrpoints = 0,
              ylim = c(-5,20), xlim = c(0,1),
              xlab = "expression specificity",
              ylab = "max expression across facets ")
grid (NULL,NULL, lty = 2, col = "grey60") 

gns <- c("MALAT1", "NEAT1", "XIST", "HAR1A", "BCAR4")
selGns <- which(gtexfData$HGNC_symbol %in% gns)
gnsSymb <- gtexfData$HGNC_symbol[selGns]
points(gtexEntropyTissue[selGns], gtexMaxAllLog[selGns], 
       pch = 20, col = "black")
text(gtexEntropyTissue[selGns], gtexMaxAllLog[selGns], 
     labels = gnsSymb, pos = c(4,1,4,4,4), cex = 0.8)

### Add median lines for entropy
mCodingent <- median(gtexEntropyTissue[sel], na.rm = T)
abline(v= mCodingent, col="deepskyblue3", lty=2)
mCodingexp <- median(gtexMaxAllLog[sel], na.rm = T)
abline(h= mCodingexp, col="deepskyblue3", lty=2)

par(mar=c(5, 0, 0, 2))
corners <- par("usr")
boxplot(gtexMaxAllLog[sel], ylim=c(-5,20), yaxt="n" , col="deepskyblue3" , frame=F, outline = F, boxwex = 0.3, staplewex = 0.3)
med <- median(gtexMaxAllLog[sel], na.rm = T)
text(x = corners[2]+0.15, y = med, labels = round(med, 2), col = "deepskyblue3", srt = 270)


################################################################
### Make plot for elncRNA and add genes to plot
sel <- FeatureTypes == "elncRNA"

par(mar=c(0, 5, 1.1, 0),cex = 2.2)
corners <- par("usr")
boxplot(gtexEntropyTissue[sel], horizontal=TRUE , ylim=c(0,1), xaxt="n" , col
        ="green4" , frame=F, outline = F, boxwex = 0.3, staplewex = 0.3,
        main = "e-lncRNA")
med <- median(gtexEntropyTissue[sel], na.rm = T)
text(x = med, y = corners[2]-0.3, labels = round(med,2), col = "green4", srt = 0)
par(mar=c(5, 5, 0, 0))
smoothScatter(gtexEntropyTissue[sel], gtexMaxAllLog[sel],
              colramp=colorRampPalette(c("white", "green4")),
              transformation = function(x) x,
              nrpoints = 0,
              ylim = c(-5,20), xlim = c(0,1),
              xlab = "expression specificity",
              ylab = "max expression across facets ")
grid (NULL,NULL, lty = 2, col = "grey60") 

gns <- c("EGOT", "LUCAT1", "HCCAT5", "TRERNA1", "DLX6-AS1")
selGns <- which(gtexfData$HGNC_symbol %in% gns)
gnsSymb <- gtexfData$HGNC_symbol[selGns]
points(gtexEntropyTissue[selGns], gtexMaxAllLog[selGns], 
       pch = 20, col = "black")
text(gtexEntropyTissue[selGns], gtexMaxAllLog[selGns], 
     labels = gnsSymb, pos = c(4,4,4,3), cex = 0.8)

### Add median lines for entropy
mCodingent <- median(gtexEntropyTissue[sel], na.rm = T)
abline(v= mCodingent, col="green4", lty=2)
mCodingexp <- median(gtexMaxAllLog[sel], na.rm = T)
abline(h= mCodingexp, col="green4", lty=2)

par(mar=c(5, 0, 0, 2))
corners <- par("usr")
boxplot(gtexMaxAllLog[sel], ylim=c(-5,20), yaxt="n" , col="green4" , frame=F, outline = F, boxwex = 0.3, staplewex = 0.3)
med <- median(gtexMaxAllLog[sel], na.rm = T)
text(x = corners[2]+0.15, y = med, labels = round(med, 2), col = "green4", srt = 270)


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
boxplot(gtexEntropyTissue~FeatureTypes2, col = c("red3", "purple", "cyan4", "green4"),
        main = "Specificity for 4 lncRNA groups", ylab = "expression specificity (log2, cpm")

### Add Grid
grid(col = "gray60")

### Create boxplot for Max expression measure 
boxplot(gtexMaxAllLog~FeatureTypes2, col = c("red3", "purple", "cyan4", "green4"),
        main = "Maximum expression for 4 lncRNA groups", ylab = "max expression ")

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


