### Plot markers
### Eddie Imada
###############################################################################
# load libraries
library(reshape2)
library(ggplot2)
library(SummarizedExperiment)
# load objs
load("gencodeGtex.rda")
load("~/Dropbox (MechPred)/FANTOM6/GTEX/objs/gtexEset.rda")
#extract exp matrix
genc <- assays(rse_scaled)[[1]]
cat <- exprs(gtexEset)
# check samples order
table(colData(rse_scaled)$smtsd == pData(gtexEset)$smtsd)
# Get markers expression
esr1 <- "ENSG00000091831"
krt1 <- "ENSG00000167768"
neurd1 <- "ENSG00000162992"
krt1c <- cat[rownames(cat) == krt1,]+1
krt1g <- genc[rownames(genc) == krt1,]+1
esr1c <- cat[rownames(cat) == esr1,]+1
esr1g <- genc[rownames(genc) == esr1,]+1
neuroc <- cat[rownames(cat) == neurod1,]+1
neurog <- genc[rownames(genc) == neurd1,]+1
# Get correlation between FC and GENCODE
corkrt <- cor(krt1c, krt1g)
corneuro <- cor(neuroc, neurog)
coresr <- cor(esr1c, esr1g)
# Created melted data frame
df <- cbind.data.frame(GENCODE = krt1g, 'FANTOM-CAT'= krt1c, Tissue=pData(gtexEset)$smtsd, GENE=rep("KRT1 - cor = 0.994", length(krt1g)))
df <- melt(df)
df1 <- cbind.data.frame(GENCODE = neurog, 'FANTOM-CAT'= neuroc, Tissue=pData(gtexEset)$smtsd, GENE=rep("NEUROD1 - cor = 0.999", length(neurog)))
df1 <- melt(df1)
df2 <- cbind.data.frame(GENCODE = esr1g, 'FANTOM-CAT'= esr1c, Tissue=pData(gtexEset)$smtsd, GENE=rep("ESR1 - cor = 0.999", length(neurog)))
df2 <- melt(df2)
df <- rbind.data.frame(df,df1, df2)

# Plot markers
p1 <- ggplot(df, aes(x=Tissue, y=log2(value), fill = variable)) +
    geom_boxplot(outlier.size = 0.25) +
    theme_bw() +
    theme(legend.position="left", 
          legend.title = element_blank(), 
          text = element_text(size = 14), 
          axis.text.x=element_text(angle=45, hjust=1),
          plot.margin = margin(0.5, 0.5, 0, 0, "cm")) +
    xlab("") +
    ylab("log2 expression") +
    scale_fill_manual(values = c("deeppink1", "deepskyblue1"))+
    #Spli by genes
    facet_wrap(~GENE,nrow = 3,scales = "free_y") +
    # Plot lines
    stat_summary(fun.y=median, geom="line", aes(group=variable), color = "black", alpha = 0.7)

# Save plot
ggsave(p1, height = 9, width = 13.5, dpi = 600, filename = "exprsGenes.jpg")
