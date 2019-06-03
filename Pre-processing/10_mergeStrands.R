#################################################################
### Luigi Marchionni
### Explore ReCount data based on FANTOM6 CAT annotation for GTEX


#################################################################
### Clean and set wd
rm(list=ls())
#setwd("/users/eimada/FANTOM6/")


### Load object
library(Biobase)
load("objs/fDataCAT.rda")
load("objs/TCGAPos.rda")
load("objs/TCGANeg.rda")
load("objs/uniqRangesNeg.rda")
load("objs/uniqRangesPos.rda")

fdataNeg <- fData(TCGANeg)
fdataPos <- fData(TCGAPos)
fdata <- rbind.data.frame(fdataPos, fdataNeg)

pdata <- pData(TCGANeg)

edataNeg <- exprs(TCGANeg)
edataPos <- exprs(TCGAPos)
edata <- rbind(edataPos, edataNeg)

uniqPos <- uniqRangesPos[uniqRangesPos$GeneCovFromRanges >= 0,]
uniqPos <- uniqPos[!duplicated(uniqPos$genes),]
uniqNeg <- uniqRangesNeg[uniqRangesNeg$GeneCovFromRanges >= 0,]
uniqNeg <- uniqNeg[!duplicated(uniqNeg$genes),]
uniq <- rbind.data.frame(uniqPos,uniqNeg)



miame <- new("MIAME", name="GTEX_CAT",
             title="TCGA expression data - CAT FANTOM genes - Recount2",
             lab="Marchionni-Langmead-Leek-Jaffe")


#################################################################
### Create the expressionSet
TCGAeset <- ExpressionSet(assayData <- edata,
                             phenoData = AnnotatedDataFrame(pdata),
                             featureData= AnnotatedDataFrame(fdata),
                             MIAME = miame)
annotation(TCGAeset) <- "CAT.FANTOM"


### Check validity
validObject(TCGAeset)


#################################################################
### ### Save
### save(gtexEset, file="objs/gtexEset.rda")
save(TCGAeset, file="objs/TCGAeset.rda")

q("no")
