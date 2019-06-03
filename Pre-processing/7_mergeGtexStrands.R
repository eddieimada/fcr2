#################################################################
### Luigi Marchionni
### Explore ReCount data based on FANTOM6 CAT annotation for GTEX


#################################################################
### Clean and set wd
rm(list=ls())
setwd("/users/eimada/FANTOM6/")


### Load object
library(Biobase)
load("objs/fDataCAT.rda")
load("objs/gtexEsetPos.rda")
load("objs/gtexEsetNeg.rda")
load("objs/uniqRangesNeg.rda")
load("objs/uniqRangesPos.rda")

fdataNeg <- fData(gtexEsetNeg)
fdataPos <- fData(gtexEsetPos)
fdata <- rbind.data.frame(fdataPos, fdataNeg)

# uniqPos <- uniqRangesPos[uniqRangesPos$GeneCovFromRanges < 100,]
# uniqPos <- uniqPos[!duplicated(uniqPos$genes),]
# uniqNeg <- uniqRangesNeg[uniqRangesNeg$GeneCovFromRanges < 100,]
# uniqNeg <- uniqNeg[!duplicated(uniqNeg$genes),]
# uniq <- rbind.data.frame(uniqPos,uniqNeg)
# 
# fdata <- fdata[uniq$genes,]
#################################################################

### Get expression for differrent classes of genes

### Gene classes
# par(mar=c(4,10,4,2))
# barplot(table(fdata$CAT_geneClass), horiz=TRUE, las =2, xlim = c(0,2000),
#         main = "Number of CAT gene classes with < 100 bp coverage",
#         col = "cornflowerblue")
# par(mar=c(4,8,4,2))
# barplot(table(fdata$CAT_DHS_type), horiz=TRUE, las =2, col = "cornflowerblue",
#         main = "Number of CAT DHS type with < 100 bp coverage")
# 
# table(fdata$CAT_geneClass,
#       fdata$CAT_DHS_type)


pdata <- pData(gtexEsetNeg)

edataNeg <- exprs(gtexEsetNeg)
edataPos <- exprs(gtexEsetPos)
edata <- rbind(edataPos, edataNeg)

uniqPos <- uniqRangesPos[uniqRangesPos$GeneCovFromRanges >= 0,]
uniqPos <- uniqPos[!duplicated(uniqPos$genes),]
uniqNeg <- uniqRangesNeg[uniqRangesNeg$GeneCovFromRanges >= 0,]
uniqNeg <- uniqNeg[!duplicated(uniqNeg$genes),]
uniq <- rbind.data.frame(uniqPos,uniqNeg)



miame <- new("MIAME", name="GTEX_CAT",
             title="GTEX expression data - CAT FANTOM genes - Recount2",
             lab="Marchionni-Langmead-Leek-Jaffe")


#################################################################
### Create the expressionSet
gtexEset <- ExpressionSet(assayData <- edata,
                             phenoData = AnnotatedDataFrame(pdata),
                             featureData= AnnotatedDataFrame(fdata),
                             MIAME = miame)
annotation(gtexEsetNeg) <- "CAT.FANTOM"


### Check validity
validObject(gtexEset)


#################################################################
### ### Save
### save(gtexEset, file="objs/gtexEset.rda")
save(gtexEset, file="objs/gtexEset.rda")


rownames(uniq) <- uniq$genes
uniq <- uniq[rownames(edata),]

newdat <- apply(edata, 2, function(x,y){
        RPK <- x/(y$GeneCovFromRanges/1000)
}, y = uniq)

png(file="figs/expressionBycov.png", height = 1000, width = 2000)
plot(uniq$GeneCovFromRanges, newdat[,1], xlim=c(0,10000), ylim=c(0,100), ylab= "normalized expression", xlab="bases covered")
dev.off()


q("no")
