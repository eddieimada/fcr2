###########################################################################
### Count DEG from TCGA 
###########################################################################

### Libraries
library(reshape2)
library(xtable)

### Set working directory
setwd("~/DropboxMech/FANTOM6/TCGA/")

### Clean environment 
rm(list=ls())

### Load
files <- list.files("~/Dropbox (MechPred)/FANTOM6/TCGA/DGE/", recursive = TRUE,
                    pattern="DGE.+rda$", full.names = TRUE)

### Remove unwanted files
files <- files[-grep("_sep", files)]
files <- files[-grep("Rectum", files)]
files <- files[-grep("Gleason", files)]
files <- files[-grep("Colon", files)]

### Read files
DGE <- lapply(files, function(x) mget(load(x))) 
DGE <- unlist(DGE, recursive = FALSE)

### Add names
nmsT <- gsub("__", "_", 
             gsub("\\.rda", "", gsub("/objs/DGE", "_",
                                     gsub(".+TCGA/DGE//", "", files))))
nmsT <- gsub("list", "", nmsT)

### Tissue names
nmsT <- gsub("CodingmRNA", "_mrna", nmsT)
nmsT <- gsub("Coding", "_mrna", nmsT)
nmsT <- gsub("mrnamRNA", "_mrna", nmsT)
nmsT <- gsub("elncRNA", "_enhc", nmsT)
nmsT <- gsub("intergenicP", "_iProm", nmsT)
nmsT <- gsub("divergentP", "_dProm", nmsT)
nmsT <- gsub("Otherlnc", "_Otherlnc", nmsT)
nmsT <- gsub("Pseudo", "_Pseudo", nmsT)
nmsT <- gsub("RNA", "_RNA", nmsT)
nmsT <- gsub("Sense_overlap", "_Senseoverlap", nmsT)

### Gene class names
nmsG <- sapply(nmsT, function(x) {
    out <- unique(unlist(strsplit(split="_", x)))
    out[ out  %in% c("mrna", "enhc", "iProm", "dProm", "Otherlnc",
                     "Pseudo", "RNA", "Senseoverlap")]
})

### Tissue
nmsT <- gsub("_.+", "", nmsT)

### Final names
nms <- paste(nmsT, nmsG, sep="_")
names(DGE) <- nms


###########################################################################
###########################################################################
### Read Gene annotaiton
ann <- read.table("~/DropboxMech/FANTOM6/GRangesAnnotation/data/F6_CAT.gene.info.tsv.gz", sep="\t",
		  header=TRUE, stringsAsFactors=FALSE, fill=TRUE, comment.char="#", quote="")

### Make table
colSel <- c("geneID", "geneName", "geneType", "CAT_geneClass", "CAT_DHS_type", "cntg", "geneStart", "geneEnd", "strnd")
### Make table
names(colSel)<- c("geneID", "geneSymbol", "geneType", "CAT_geneClass", "CAT_DHS_type",
		  "Chromosome", "geneStart", "geneEnd", "strand")


### To print long tables prepare elements for longtable
addtorow <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(0)
addtorow$command  <- c(paste("\\hline \n", "\\endhead \n", "\\hline \n",
                             "{\\tiny Continued on next page} \n", "\\endfoot \n",
			     "\\endlastfoot \n", sep=""))


##########################################################################
##########################################################################
###  coding mRNA
##########################################################################

### Subset mrna
cancerType <- grep("mrna", names(DGE), value=TRUE)
tmp <- DGE[cancerType]

### Add cancer type
tmp <- mapply(dat=tmp, nms=names(tmp), function(dat, nms) {
	dat$cancerType <- nms
	dat
},  SIMPLIFY = FALSE)

### Merge to correct the p.value globally
tmp <- do.call("rbind",  tmp)

### Check dimensions
dim(tmp)

### Correct for multiple testing: BH
tmp$globalBH <- p.adjust(tmp$P.Value, method="BH")
### Correct for multiple testing: Bonferroni
tmp$globalBONF <- p.adjust(tmp$P.Value, method="bonferroni")

### Select UP-regulate only 
selUP <- tmp$t > 0
### UP-regulate only using BH
tmp$globalBHup <- NA
tmp$globalBHup[selUP] <- p.adjust(tmp$P.Value[selUP], method="BH")
### UP-regulate only  using Bonferroni
tmp$globalBONFup <- NA
tmp$globalBONFup[selUP] <- p.adjust(tmp$P.Value[selUP], method="bonferroni")

### Select DOWN-regulate only 
selDW <- tmp$t < 0
### DW-regulate only using BH
tmp$globalBHdw <- NA
tmp$globalBHdw[selDW] <- p.adjust(tmp$P.Value[selDW], method="BH")
### DW-regulate only using Bonferroni
tmp$globalBONFdw <- NA
tmp$globalBONFdw[selDW] <- p.adjust(tmp$P.Value[selDW], method="bonferroni")

### ### Count BH
### table(tmp$globalBH < 0.01)
### table(SIG=tmp$globalBH < 0.01, CANCER=tmp$cancerType)
### table(SIG=tmp$globalBHup < 0.01, CANCER=tmp$cancerType)

### Split back
tmp <- lapply(cancerType, function(nms,  dat) {
	dat[dat$cancerType %in% nms,  ]
},  dat=tmp)
names(tmp) <- cancerType

### Check
sapply(tmp, dim)
tmp <- lapply(tmp, function(x){
    x[order(abs(x$logFC), decreasing = T),]
})

### Intersection for UP-regulated with BH
upGenesBH <- lapply(tmp, function(x) x$geneName[ x$t > 0 & x$globalBHup < 0.000001  & !is.na(x$globalBHup) ] )
upGenesBHmrna <- Reduce("intersect", upGenesBH)
str(upGenesBHmrna)

### Intersection for DW-regulated with BH
dwGenesBH <- lapply(tmp, function(x) x$geneName[ x$t < 0 & x$globalBHdw < 0.000001  & !is.na(x$globalBHdw) ] )
dwGenesBHmrna <- Reduce("intersect", dwGenesBH)
str(dwGenesBHmrna)

### Intersection for UP-regulated with Bonferroni
upGenesBONF <- lapply(tmp, function(x) x$geneName[ x$t > 0 & x$globalBONFup < 0.01  & !is.na(x$globalBONFup) ] )
upGenesBONFmrna <- Reduce("intersect", upGenesBONF)
str(upGenesBONFmrna)

### Intersection for DW-regulated with Bonferroni
dwGenesBONF <- lapply(tmp, function(x) x$geneName[ x$t < 0 & x$globalBONFdw < 0.01  & !is.na(x$globalBONFdw) ] )
dwGenesBONFmrna <- Reduce("intersect", dwGenesBONF)
str(dwGenesBONFmrna)


### Make table
tmpUp <- ann[ann$geneName %in% upGenesBHmrna, colSel]
colnames(tmpUp) <- names(colSel)
tmpUp$Direction <- "Up-regulated in tumor"
tmpDw <- ann[ann$geneName %in% dwGenesBHmrna, colSel]
colnames(tmpDw) <- names(colSel)
tmpDw$Direction <- "Down-regulated in tumor"
tmp <- rbind(tmpUp, tmpDw)

### Write table
write.csv(tmp, file="./textFinal/mRNA_geneIntersection.csv", row.names=FALSE)

### Write a latexTable
sink(file="./textFinal/mRNA_geneIntersection.tex")
print(xtable(tmp,
	     caption = "Intersection of significant mRNA across 13 cancer types (Global FDR < 0.000001)",
	     label = "tab:mRNAintersection"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=FALSE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)




##########################################################################
##########################################################################
### dProm
##########################################################################

### Subset dProm
cancerType <- grep("dProm", names(DGE), value=TRUE)
tmp <- DGE[cancerType]

### Add cancer type
tmp <- mapply(dat=tmp, nms=names(tmp), function(dat, nms) {
	dat$cancerType <- nms
	dat
},  SIMPLIFY = FALSE)

### Merge to correct the p.value globally
tmp <- do.call("rbind",  tmp)

### Check dimensions
dim(tmp)

### Correct for multiple testing: BH
tmp$globalBH <- p.adjust(tmp$P.Value, method="BH")
### Correct for multiple testing: Bonferroni
tmp$globalBONF <- p.adjust(tmp$P.Value, method="bonferroni")

### Select UP-regulate only 
selUP <- tmp$t > 0
### UP-regulate only using BH
tmp$globalBHup <- NA
tmp$globalBHup[selUP] <- p.adjust(tmp$P.Value[selUP], method="BH")
### UP-regulate only  using Bonferroni
tmp$globalBONFup <- NA
tmp$globalBONFup[selUP] <- p.adjust(tmp$P.Value[selUP], method="bonferroni")

### Select DOWN-regulate only 
selDW <- tmp$t < 0
### DW-regulate only using BH
tmp$globalBHdw <- NA
tmp$globalBHdw[selDW] <- p.adjust(tmp$P.Value[selDW], method="BH")
### DW-regulate only using Bonferroni
tmp$globalBONFdw <- NA
tmp$globalBONFdw[selDW] <- p.adjust(tmp$P.Value[selDW], method="bonferroni")

### ### Count BH
### table(tmp$globalBH < 0.01)
### table(SIG=tmp$globalBH < 0.01, CANCER=tmp$cancerType)
### table(SIG=tmp$globalBHup < 0.01, CANCER=tmp$cancerType)

### Split back
tmp <- lapply(cancerType, function(nms,  dat) {
	dat[dat$cancerType %in% nms,  ]
},  dat=tmp)
names(tmp) <- cancerType

### Check
sapply(tmp, dim)
tmp <- lapply(tmp, function(x){
    x[order(abs(x$logFC), decreasing = T),]
})

### Intersection for UP-regulated with BH
upGenesBH <- lapply(tmp, function(x) x$geneName[ x$t > 0 & x$globalBHup < 0.1  & !is.na(x$globalBHup) ] )
upGenesBHdProm <- Reduce("intersect", upGenesBH)
str(upGenesBHdProm)

### Intersection for DW-regulated with BH
dwGenesBH <- lapply(tmp, function(x) x$geneName[ x$t < 0 & x$globalBHdw < 0.1  & !is.na(x$globalBHdw) ] )
dwGenesBHdProm <- Reduce("intersect", dwGenesBH)
str(dwGenesBHdProm)

### Intersection for UP-regulated with Bonferroni
upGenesBONF <- lapply(tmp, function(x) x$geneName[ x$t > 0 & x$globalBONFup < 0.1  & !is.na(x$globalBONFup) ] )
upGenesdBONFProm <- Reduce("intersect", upGenesBONF)
str(upGenesdBONFProm)

### Intersection for DW-regulated with Bonferroni
dwGenesBONF <- lapply(tmp, function(x) x$geneName[ x$t < 0 & x$globalBONFdw < 0.1  & !is.na(x$globalBONFdw) ] )
dwGenesdBONFProm <- Reduce("intersect", dwGenesBONF)
str(dwGenesdBONFProm)

### Make table
tmpUp <- ann[ann$geneName %in% upGenesBHdProm, colSel]
colnames(tmpUp) <- names(colSel)
tmpUp$Direction <- "Up-regulated in tumor"
tmpDw <- ann[ann$geneName %in% dwGenesBHdProm, colSel]
colnames(tmpDw) <- names(colSel)
tmpDw$Direction <- "Down-regulated in tumor"
tmp <- rbind(tmpUp, tmpDw)

### Write table
write.csv(tmp, file="./textFinal/dProm_geneIntersection.csv", row.names=FALSE)

### Write a latexTable
sink(file="./textFinal/dProm_geneIntersection.tex")
print(xtable(tmp,
	     caption = "Intersection of significant divergent promoters across 13 cancer types (Global FDR < 0.1)",
	     label = "tab:dpromRNAintersection"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=FALSE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)



#########################################################################
#########################################################################
### iProm
##########################################################################

### Subset iProm
cancerType <- grep("iProm", names(DGE), value=TRUE)
tmp <- DGE[cancerType]

### Add cancer type
tmp <- mapply(dat=tmp, nms=names(tmp), function(dat, nms) {
	dat$cancerType <- nms
	dat
},  SIMPLIFY = FALSE)

### Merge to correct the p.value globally
tmp <- do.call("rbind",  tmp)

### Check dimensions
dim(tmp)

### Correct for multiple testing: BH
tmp$globalBH <- p.adjust(tmp$P.Value, method="BH")
### Correct for multiple testing: Bonferroni
tmp$globalBONF <- p.adjust(tmp$P.Value, method="bonferroni")

### Select UP-regulate only 
selUP <- tmp$t > 0
### UP-regulate only using BH
tmp$globalBHup <- NA
tmp$globalBHup[selUP] <- p.adjust(tmp$P.Value[selUP], method="BH")
### UP-regulate only  using Bonferroni
tmp$globalBONFup <- NA
tmp$globalBONFup[selUP] <- p.adjust(tmp$P.Value[selUP], method="bonferroni")

### Select DOWN-regulate only 
selDW <- tmp$t < 0
### DW-regulate only using BH
tmp$globalBHdw <- NA
tmp$globalBHdw[selDW] <- p.adjust(tmp$P.Value[selDW], method="BH")
### DW-regulate only using Bonferroni
tmp$globalBONFdw <- NA
tmp$globalBONFdw[selDW] <- p.adjust(tmp$P.Value[selDW], method="bonferroni")

### ### Count BH
### table(tmp$globalBH < 0.01)
### table(SIG=tmp$globalBH < 0.01, CANCER=tmp$cancerType)
### table(SIG=tmp$globalBHup < 0.01, CANCER=tmp$cancerType)

### Split back
tmp <- lapply(cancerType, function(nms,  dat) {
	dat[dat$cancerType %in% nms,  ]
},  dat=tmp)
names(tmp) <- cancerType

### Check
sapply(tmp, dim)
tmp <- lapply(tmp, function(x){
    x[order(abs(x$logFC), decreasing = T),]
})

### Intersection for UP-regulated with BH
upGenesBH <- lapply(tmp, function(x) x$geneName[ x$t > 0 & x$globalBHup < 0.1  & !is.na(x$globalBHup) ] )
upGenesBHiProm <- Reduce("intersect", upGenesBH)
str(upGenesBHiProm)

### Intersection for DW-regulated with BH
dwGenesBH <- lapply(tmp, function(x) x$geneName[ x$t < 0 & x$globalBHdw < 0.1  & !is.na(x$globalBHdw) ] )
dwGenesBHiProm <- Reduce("intersect", dwGenesBH)
str(dwGenesBHiProm)

### Intersection for UP-regulated with Bonferroni
upGenesBONF <- lapply(tmp, function(x) x$geneName[ x$t > 0 & x$globalBONFup < 0.1  & !is.na(x$globalBONFup) ] )
upGenesBONFiProm <- Reduce("intersect", upGenesBONF)
str(upGenesBONFiProm)

### Intersection for DW-regulated with Bonferroni
dwGenesBONF <- lapply(tmp, function(x) x$geneName[ x$t < 0 & x$globalBONFdw < 0.1  & !is.na(x$globalBONFdw) ] )
dwGenesBONFiProm <- Reduce("intersect", dwGenesBONF)
str(dwGenesBONFiProm)

### Make table
tmpUp <- ann[ann$geneName %in% upGenesBHiProm, colSel]
colnames(tmpUp) <- names(colSel)
tmpUp$Direction <- "Up-regulated in tumor"
tmpDw <- ann[ann$geneName %in% dwGenesBHiProm, colSel]
colnames(tmpDw) <- names(colSel)
tmpDw$Direction <- "Down-regulated in tumor"
tmp <- rbind(tmpUp, tmpDw)

### Write table
write.csv(tmp, file="./textFinal/iProm_geneIntersection.csv", row.names=FALSE)

### Write a latexTable
sink(file="./textFinal/iProm_geneIntersection.tex")
print(xtable(tmp,
             caption = "Intersection of significant intergenic promoters across 13 cancer types (Global FDR < 0.1)",
             label = "tab:ipromRNAintersection"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=FALSE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)


##########################################################################
##########################################################################
### enhc
##########################################################################

### Subset enhancers
cancerType <- grep("enhc", names(DGE), value=TRUE)
tmp <- DGE[cancerType]

### Add cancer type
tmp <- mapply(dat=tmp, nms=names(tmp), function(dat, nms) {
	dat$cancerType <- nms
	dat
},  SIMPLIFY = FALSE)

### Merge to correct the p.value globally
tmp <- do.call("rbind",  tmp)

### Check dimensions
dim(tmp)

### Correct for multiple testing: BH
tmp$globalBH <- p.adjust(tmp$P.Value, method="BH")
### Correct for multiple testing: Bonferroni
tmp$globalBONF <- p.adjust(tmp$P.Value, method="bonferroni")

### Select UP-regulate only 
selUP <- tmp$t > 0
### UP-regulate only using BH
tmp$globalBHup <- NA
tmp$globalBHup[selUP] <- p.adjust(tmp$P.Value[selUP], method="BH")
### UP-regulate only  using Bonferroni
tmp$globalBONFup <- NA
tmp$globalBONFup[selUP] <- p.adjust(tmp$P.Value[selUP], method="bonferroni")

### Select DOWN-regulate only 
selDW <- tmp$t < 0
### DW-regulate only using BH
tmp$globalBHdw <- NA
tmp$globalBHdw[selDW] <- p.adjust(tmp$P.Value[selDW], method="BH")
### DW-regulate only using Bonferroni
tmp$globalBONFdw <- NA
tmp$globalBONFdw[selDW] <- p.adjust(tmp$P.Value[selDW], method="bonferroni")

### ### Count BH
### table(tmp$globalBH < 0.01)
### table(SIG=tmp$globalBH < 0.01, CANCER=tmp$cancerType)
### table(SIG=tmp$globalBHup < 0.01, CANCER=tmp$cancerType)

### Split back
tmp <- lapply(cancerType, function(nms,  dat) {
	dat[dat$cancerType %in% nms,  ]
},  dat=tmp)
names(tmp) <- cancerType

### Check
sapply(tmp, dim)
tmp <- lapply(tmp, function(x){
    x[order(abs(x$logFC), decreasing = T),]
})
### Intersection for UP-regulated with BH
upGenesBH <- lapply(tmp, function(x) x$geneName[ x$t > 0 & x$globalBHup < 0.1  & !is.na(x$globalBHup) ] )
upGenesBHelnc <- Reduce("intersect", upGenesBH)
str(upGenesBHelnc)

### Intersection for DW-regulated with BH
dwGenesBH <- lapply(tmp, function(x) x$geneName[ x$t < 0 & x$globalBHdw < 0.1  & !is.na(x$globalBHdw) ] )
dwGenesBHelnc <- Reduce("intersect", dwGenesBH)
str(dwGenesBHelnc)

### Intersection for UP-regulated with Bonferroni
upGenesBONF <- lapply(tmp, function(x) x$geneName[ x$t > 0 & x$globalBONFup < 0.1  & !is.na(x$globalBONFup) ] )
upGenesBONFelnc <- Reduce("intersect", upGenesBONF)
str(upGenesBONFelnc)

### Intersection for DW-regulated with Bonferroni
dwGenesBONF <- lapply(tmp, function(x) x$geneName[ x$t < 0 & x$globalBONFdw < 0.1  & !is.na(x$globalBONFdw) ] )
dwGenesBONFelnc <- Reduce("intersect", dwGenesBONF)
str(dwGenesBONFelnc)


### Make table
tmpUp <- ann[ann$geneName %in% upGenesBHelnc, colSel]
colnames(tmpUp) <- names(colSel)
tmpUp$Direction <- "Up-regulated in tumor"
tmpDw <- ann[ann$geneName %in% dwGenesBHelnc, colSel]
colnames(tmpDw) <- names(colSel)
tmpDw$Direction <- "Down-regulated in tumor"
tmp <- rbind(tmpUp, tmpDw)

### Write table
write.csv(tmp, file="./textFinal/elncRNA_geneIntersection.csv", row.names=FALSE)

### Write a latexTable
sink(file="./textFinal/elncRNA_geneIntersection.tex")
print(xtable(tmp,
             caption = "Intersection of significant enhancers across 13 cancer types (Global FDR < 0.1)",
             label = "tab:elncRNAintersection"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=FALSE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

###########################################################################
### Saving up- and down-regulated genes by RNA subtypes.

BHlistUP <- list(dProm=upGenesBHdProm,elnc=upGenesBHelnc, iProm=upGenesBHiProm, mrna=upGenesBHmrna)
BHlistDN <- list(dProm=dwGenesBHdProm,elnc=dwGenesBHelnc, iProm=dwGenesBHiProm, mrna=dwGenesBHmrna)

save(BHlistUP, BHlistDN, file = "./objs/BHintersectList.rda")
###########################################################################
### Session information and clean quit

### Session
sessionInfo()

### Clean
rm(list=ls())

### Quit
q("no")
