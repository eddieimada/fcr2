#################################################################
###Tejasvi Matam
### Create plot of % genes expressed/facet for each lncRNA type

#################################################################
### Clean and set wd
rm(list=ls())
#################################################################
### Load packages & data
require(Biobase)
require(plotrix)
require(edgeR)
# Path to GTEx object
path <- "~/FANTOM6/GTEX/gtexEset_raw.rda"
load(path)
load("~/FANTOM6/GTEX/uniqRangesNeg.rda")
load("~/FANTOM6/GTEX/uniqRangesPos.rda")
#################################################################
### Get gene length
uniqPos <- uniqRangesPos[!duplicated(uniqRangesPos$genes),]
uniqNeg <- uniqRangesNeg[!duplicated(uniqRangesNeg$genes),]
uniq <- rbind.data.frame(uniqPos,uniqNeg)
rownames(uniq) <- uniq$genes

### Create TPM matrix
gtexExprs<- exprs(gtexEset_raw)
uniq <- uniq[rownames(gtexExprs),]
KM <- uniq$GeneCovFromRanges

gtexTPM <- sweep(gtexTPM, 1, FUN = "/", STATS = KM)

libsize <- colSums(gtexTPM)
normLib <- libsize/1000000

gtexTPM <- sweep(gtexTPM, 2, FUN = "/", STATS =  normLib)

#################################################################
### New matrix of 0's and 1's based on TPM 
gtexNewexprs<- gtexTPM >= 0.01

gtexCount<- gtexNewexprs*1 

#################################################################
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
 
##################################################################
### Sum types by FeatureTypes
gtexCountSum <- (apply(gtexCount, 2, function(x, y) 
    {tapply(x,y, sum, na.rm= TRUE)}, y= FeatureTypes)) 

### Compute percent of genes expressed from total genes
if (all(rownames(gtexCountSum) == names(table(FeatureTypes))) ){
   gtexCountpercent<- apply(gtexCountSum, 2, function(x, y) {
     (x/y)*100
     }, y = table(FeatureTypes))
   print("All is good")
 } else print("Check order")

#################################################################
### Collapse by tissue types 
TissueTypes <- factor(gtexEset$smtsd)

#################################################################
### Find quantiles using TissueType
gtexSummary <- (apply(gtexCountpercent, 1, function(x, y) {
  ci05 <- tapply(x,y, quantile, c(0.05), na.rm= TRUE)
  ci95 <- tapply(x,y, quantile, c(0.95), na.rm= TRUE)
  m <- tapply(x,y, mean, na.rm= TRUE)
  out <- list(CI05=ci05,m=m,CI95=ci95)
  }, y=TissueTypes)) 

#################################################################
### Rearrange from highest to lowest
newOrder2 <- order(gtexSummary$elncRNA$m, decreasing = TRUE)

### Apply new order to all 
gtexMean2 <- gtexSummary$elncRNA$m[newOrder2]
gtexCI05order2 <- gtexSummary$elncRNA$CI05[newOrder2]
gtexCI95order2 <- gtexSummary$elncRNA$CI95[newOrder2]

gtexMean3 <- gtexSummary$`divergent p-lncRNA`$m[newOrder2]
gtexCI05order3 <- gtexSummary$`divergent p-lncRNA`$CI05[newOrder2]
gtexCI95order3 <- gtexSummary$`divergent p-lncRNA`$CI95[newOrder2]

gtexMean4 <- gtexSummary$`intergenic p-lncRNA`$m[newOrder2]
gtexCI05order4 <- gtexSummary$`intergenic p-lncRNA`$CI05[newOrder2]
gtexCI95order4 <- gtexSummary$`intergenic p-lncRNA`$CI95[newOrder2]

gtexMean5 <- gtexSummary$codingmRNA$m[newOrder2]
gtexCI05order5 <- gtexSummary$codingmRNA$CI05[newOrder2]
gtexCI95order5 <- gtexSummary$codingmRNA$CI95[newOrder2]

#################################################################
###Save objects

save(gtexTPM, file = "./objs/gtexTPM.rda")

save(gtexCountSum, file = "./objs/gtexCountSum.rda")

save(gtexCount, file = "./objs/gtexCount.rda")

save(gtexSummary, file = "./objs/gtexSummary.rda")

#################################################################
### Save Image 
png(filename = "./figs/percent_genes_expressed_by_tissue_type.png", 
      width=6000, height=3000, res=300)

par(mar = c(9,4,1,1))

plotCI(x = gtexMean2,
       ui = gtexCI95order2, sfrac = 0.005,
       li = gtexCI05order2, scol = "grey70",
       ylab = "mean % of genes expressed in each facet", xlab = "", xaxt = "n",
       col = "green4", ylim = c(0,100),
       pch=16)

mEln <- mean(gtexSummary$elncRNA$m)
abline(h=mEln, col="green4", lty=2)
text(55.1, mEln+2, round(mEln, 2), col="green4", cex = .7)


plotCI(x = gtexMean3,
       ui = gtexCI95order3, sfrac = 0.005,
       li = gtexCI05order3, scol = "grey70",
       ylab = "mean % of genes expressed in each facet", xlab = "", xaxt = "n",
       col = "purple3", ylim = c(0,100),
       pch=16,add = TRUE)

mDivP <- mean(gtexSummary$`divergent p-lncRNA`$m)
abline(h=mDivP, col="purple3", lty=2)
text(55.1, mDivP+2, round(mDivP, 2), col="purple3", cex = .7)


plotCI(x = gtexMean4,
       ui = gtexCI95order4, sfrac = 0.005,
       li = gtexCI05order4, scol = "grey70",
       ylab = "mean % of genes expressed in each facet", xlab = "", xaxt = "n",
       col = "deepskyblue3", ylim = c(0,100),
       pch=16,add = TRUE)

mInter <- mean(gtexSummary$`intergenic p-lncRNA`$m)
abline(h=mInter, col="deepskyblue3", lty=2)
text(55.1, mInter+2, round(mInter, 2), col="deepskyblue3", cex = .7)


plotCI(x = gtexMean5,
       ui = gtexCI95order5, sfrac = 0.005,
       li = gtexCI05order5, scol = "grey70",
       ylab = "mean % of genes expressed in each facet", xlab = "", xaxt = "n",
       col = "firebrick1", ylim = c(0,100),
       pch=16,add = TRUE)

mCoding <- mean(gtexSummary$codingmRNA$m)
abline(h=mCoding, col="firebrick1", lty=2)
text(55.1,mCoding+2, round(mCoding, 2), col="firebrick1", cex = .7)

### Add x-axis facet groups 
staxlab(1, at = 1:54, labels = rownames(gtexMean2), cex.axis = .5, las = 3, srt = 45)

### Add grid 
grid(col= "grey60")

### Create legend 
legend("bottomleft",  
       c("elncRNA", "divergent p-lncRNA", "intergenic p-lncRNA", "codingmRNA"), 
       fill = c("green4", "purple3", "deepskyblue3", "firebrick1"), cex = 1, bty = "n")

dev.off()

#################################################################
### Session information and clean, quit

### Session
sessionInfo()

### Clean
rm(list=ls())

### Quit
q("no")


