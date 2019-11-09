### Load libraries & objs
library(Biobase)
library(ggplot2)
load("/Users/elimada/Dropbox (MechPred)/FANTOM6/GTEX/objs/gtexTPM.rda")
load("/Users/elimada/Dropbox (MechPred)/FANTOM6/GTEX/objs/gtexEset.rda")
### extract matrices
mat <- gtexCPM
rm(gtexCPM)
gc()
feat <- fData(gtexEset)
pheno <- pData(gtexEset)
rm(gtexEset)
gc()

### get counts of expressed genes
fc.sum <- apply(mat, 1, function(x){
    sum(x > 0.01)
})

### get percentage of samples expressed in gtex
fc.pc <- fc.sum/9662

### Plot histograms
pcs <- qplot(fc.pc, geom = "histogram", bins = 200) + 
    theme_bw() +
    xlab("Percentage of samples")
ggsave(pcs, height = 6, width = 8, dpi = 300, filename = "~/Dropbox (MechPred)/FANTOM6/Manuscripts/FCR2/review/pcSamples.png")
pcs.zoom <- qplot(fc.pc[fc.pc < 0.01], bins=200) +
    theme_bw() +
    xlab("Percentage of samples")
ggsave(pcs.zoom, height = 6, width = 8, dpi = 300, filename = "~/Dropbox (MechPred)/FANTOM6/Manuscripts/FCR2/review/pcSamples_zoom.png")

### get genes not expressed & plot distribution by class
gns <- names(fc.pc[fc.pc < 0.01])
length(gns)
table(grepl("^CAT", gns))
perm <- table(names(fc.pc) %in% gns)
gns.class <- table(feat[gns,]$RNAclass)
gns.class <- gns.class[-5]
gns.class.pc <- gns.class/length(gns)
df <- data.frame(gns.class.pc)
df <- df[order(df$Freq, decreasing = T),]
df$Var1 <- factor(df$Var1, df$Var1)
plot <- ggplot(df, aes(x = Var1, y = Freq)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("RNA class") +
    ylab("Percentage")

ggsave(plot, height = 6, width = 8, dpi = 300, filename = "~/Dropbox (MechPred)/FANTOM6/Manuscripts/FCR2/review/pcRNAclass.png")


### perform chi-square tests according to different FANTOM-CAT sets. Files bellow can be obtained from the FANTOM-CAT website.
robust <- read.csv("~/Downloads/FANTOM_CAT.lv3_robust.bed.gz", sep = "\t", header = F)
rgns <- unique(gsub("\\..+", "",robust$V4))
rne <- rgns[rgns %in% gns]
rob <- table(rgns %in% rne)

stringent <- read.csv("~/Downloads/FANTOM_CAT.lv4_stringent.bed.gz", sep = "\t", header = F)
sgns <- unique(gsub("\\..+", "",stringent$V4))
sne <- sgns[sgns %in% gns]
stringt <- table(sgns %in% sne)

chi <- cbind(perm,rob, stringt)

chisq.test(chi)
