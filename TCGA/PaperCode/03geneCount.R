###########################################################################
### Count DEG from TCGA 
###########################################################################

### Libraries
library(reshape2)

### Set working directory
setwd("~/Research/FANTOM6/ReCount2/catTCGA/")

### Clean environment 
rm(list=ls())

### Load
files <- list.files("~/DropboxMech/FANTOM6/TCGA/DGE/", recursive = TRUE,
                    pattern="DGE.+csv$", full.names = TRUE)
files <- files[-grep("Unscaled", files)]

### Read files
DGE <- lapply(files, read.csv, stringsAsFactors=FALSE)

### Add names
nmsT <- gsub("__", "_", 
                   gsub("\\.csv", "", gsub("/text/DGE", "_",
                                           gsub(".+TCGA/DGE//", "", files))))

### Tissue names
nmsT <- gsub("CodingmRNA", "_mrna", nmsT)
nmsT <- gsub("Coding", "_mrna", nmsT)
nmsT <- gsub("mrnamRNA", "_mrna", nmsT)
nmsT <- gsub("elncRNA", "_enhc", nmsT)
nmsT <- gsub("intergenicP", "_iProm", nmsT)
nmsT <- gsub("divergentP", "_dProm", nmsT)

### Gene class names
nmsG <- sapply(nmsT, function(x) {
  out <- unique(unlist(strsplit(split="_", x)))
  out[ out  %in% c("mrna", "enhc", "iProm", "dProm")]
})

### Tissue
nmsT <- gsub("_.+", "", nmsT)

### Final names
nms <- paste(nmsT, nmsG, sep="_")
names(DGE) <- nms

###########################################################################
### filter
countDGE <- t(sapply(DGE, function(x) {
  if (all(is.na(x$t))) c("FALSE"=0, "TRUE"=0)
  else table(x$t > 0 & x$adj.P.Val < 0.01)
  }))

### Groups by cancer type
tissue <- gsub("_.+", "", rownames(countDGE))
table(tissue)

### UP dge
datUP <- matrix(countDGE[, "TRUE"], ncol = 4, byrow = T)
colnames(datUP) <- c("coding", "divP", "enhc", "intP")
rownames(datUP) <- unique(tissue)
datUP
write.csv(datUP, file="text/datUPconsensus.csv")

### DW dge
datDW <- matrix(countDGE[, "FALSE"], ncol = 4, byrow = T)
colnames(datDW) <- c("coding", "divP", "enhc", "intP")
rownames(datDW) <- unique(tissue)
datDW
write.csv(datDW, file="text/datDWconsensus.csv")


#################################################################
### Up genes
upGns <- lapply(DGE, function(x){
  sel <- x$t > 0 & x$adj.P.Val < 0.01
  paste(x$geneID[sel], gsub("\\..$", "", x$geneName[sel]), sep  ="_")
})

Reduce("intersect", upGns[grep("_mrna", names(upGns))])
Reduce("intersect", upGns[grep("_enhc", names(upGns))])
Reduce("intersect", upGns[grep("_dProm", names(upGns))])
Reduce("intersect", upGns[grep("_iProm", names(upGns))])


#################################################################
### dw genes
dwGns <- lapply(DGE, function(x){
  sel <- x$t < 0 & x$adj.P.Val < 0.01
  paste(x$geneID[sel], gsub("\\..$", "", x$geneName[sel]), sep  ="_")
})

Reduce("intersect", dwGns[grep("_mrna", names(dwGns))])
Reduce("intersect", dwGns[grep("_enhc", names(dwGns))])
Reduce("intersect", dwGns[grep("_dProm", names(dwGns))])
Reduce("intersect", dwGns[grep("_iProm", names(dwGns))])


#################################################################
### Session information and clean quit

### Session
sessionInfo()

### Clean
rm(list=ls())

### Quit
q("no")