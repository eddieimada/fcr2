### Plot markers
### Eddie Imada
###############################################################################
# Load libraries
library(reshape2)
library(ggplot2)
library(SummarizedExperiment)
# Set wd
setwd("~/Dropbox (MechPred)/FANTOM6/GTEX/")
# select markers
sel <- c("CD4", "MYH11", "ACTA2", "INS", "INSR", "KRT2", "KRT3")

# load dataests
load("~/gencodeGtex.rda")
load("~/Dropbox (MechPred)/FANTOM6/GTEX/objs/gtexEset.rda")
# extract exp matrices and other infos
cat <- exprs(gtexEset)
genc <- assays(rse_scaled)[[1]]
pheno <- pData(gtexEset)
feat <- fData(gtexEset)
# Remove objs to clear memory
rm(rse_scaled, gtexEset)
gc()

# Crate GENE_ID to SYMBOL mapping vector
symbol <- setNames(feat$geneName, feat$geneID)

# Select markers
keep <- match(sel, feat$geneName)
keep2 <- match(feat$geneID[keep], rownames(genc))
cat <- cat[keep,]
genc <-genc[keep2,]
# Rename rows to symbols
rownames(cat) <- symbol[rownames(cat)]
rownames(genc) <- symbol[rownames(genc)]
# Create data frame to plot
df <- cbind.data.frame(t(cat)+1, Tissue=pheno$smtsd)
df2 <- cbind.data.frame(t(genc)+1, Tissue=pheno$smtsd)
df <- melt(df)
df2 <- melt(df2)
df <- cbind.data.frame(df, study=rep("FANTOM-CAT", nrow(df)))
df2 <- cbind.data.frame(df2, study=rep("GENCODE", nrow(df)))
df <- rbind.data.frame(df,df2)
# Plot markers
p1 <- ggplot(df, aes(x=Tissue, y=log2(value), fill = study)) +
    geom_boxplot(outlier.size = 0.25) +
    theme_bw() +
    theme(legend.position="left", 
          legend.title = element_blank(), 
          text = element_text(size = 14), 
          axis.text.x=element_text(angle=45, hjust=1),
          plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
          strip.text = element_text(size = 10)) +
    facet_wrap(~variable, scales = "free_y", nrow = 10, strip.position = "left") +
    xlab("") +
    ylab("log2 expression") +
    scale_fill_manual(values = c("deeppink1", "deepskyblue1")) +
    stat_summary(fun.y=median, geom="line", aes(group=study), color = "black", alpha = 0.7)

# SAve
ggsave(p1, height = 11, width = 13.5, dpi = 300, filename = "figs/extraMarkersGENCODE.jpg")
