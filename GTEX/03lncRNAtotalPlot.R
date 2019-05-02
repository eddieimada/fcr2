#################################################################
### Create plot of total number of lncRNA expressed within their tissue facet 

#################################################################
### Clean and set wd
rm(list=ls())

#################################################################
### Load packages
require(Biobase)
require(plotrix)

### Load Data 
# Path to GTEx object
path <- "~/FANTOM6/GTEX/gtexEset.rda"
load(path)
load("~/FANTOM6/GTEX/gtexCount.rda")

#################################################################
### Create lncRNA category
FeatureTypes <- factor(paste(featureData(gtexEset)$CAT_geneClass,featureData(gtexEset)$CAT_DHS_type))
levels(FeatureTypes)

### Create lncRNA group
levels(FeatureTypes)[grep("lncRNA_intergenic DHS_promoter", levels(FeatureTypes))] <- "lncRNA"
levels(FeatureTypes)[grep("lncRNA_divergent DHS_promoter", levels(FeatureTypes))] <- "lncRNA"
levels(FeatureTypes)[grep("lncRNA_divergent DHS_enhancer", levels(FeatureTypes))] <- "lncRNA" 
levels(FeatureTypes)[grep("lncRNA_intergenic DHS_enhancer", levels(FeatureTypes))] <- "lncRNA" 
levels(FeatureTypes)[grep("lncRNA_antisense DHS_enhancer", levels(FeatureTypes))] <- "lncRNA" 
levels(FeatureTypes)[grep("lncRNA_sense_intronic DHS_enhancer", levels(FeatureTypes))] <- "lncRNA" 

### Create other-lncrna category (not needed in barplot)
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

levels(FeatureTypes)

#################################################################
### Sum types by FeatureTypes
gtexCountSum <- (apply(gtexCount, 2, function(x, y) 
    {tapply(x,y, sum, na.rm= TRUE)}, y= FeatureTypes)) 

#################################################################
### Factor tissue types 
TissueTypes <- factor(gtexEset$smtsd)

#################################################################
### Collapse by tissue type
gtexSummary <- t(apply(gtexCountSum, 1, function(x, y)
  {tapply(x,y, mean, na.rm= TRUE) }, y=TissueTypes)) 


#################################################################
### Reorder from highest to lowest
newOrder <- order(gtexSummary[3,], decreasing = TRUE)
gtexSummary2 <- gtexSummary[,newOrder]


#################################################################
### Save plot 
jpeg(filename = "./figs/total numbers of lncRNA expressed.jpeg", 
     width=4000, height=1600, res=300)

par(mar = c(10,5,1,1))

barplotRNA <- barplot(gtexSummary2[3,1:54], cex.axis = .65, cex.names = .55,las = 2, ylim = c(0,25000),
ylab = "number of lncRNA genes expressed in facets", axes = TRUE)

### Add line to show mean across all samples 
lncRNAm<- mean(gtexSummary2[3,])
abline(h=lncRNAm, col="gray40", lty=2)
text(65.9, lncRNAm+900, round(lncRNAm, 2), col="black", cex = .65)

### Add grid
grid(col= "grey60")

### Get Min count 
gtexCountmin <- t(apply(gtexCountSum, 1, function(x, y)
{tapply(x,y, quantile,.05, na.rm= TRUE) }, y=TissueTypes)) 


### Get Max count 
gtexCountmax <- t(apply(gtexCountSum, 1, function(x, y)
{tapply(x,y, quantile, .95, na.rm= TRUE) }, y=TissueTypes)) 


### Plot the segments to create final barplot with bars
segments(barplotRNA, gtexCountmin[3,newOrder],
         barplotRNA, gtexCountmax[3,newOrder])

dev.off()

#################################################################
### Session information and clean, quit

### Session
sessionInfo()

### Clean
rm(list=ls())

### Quit
q("no")


