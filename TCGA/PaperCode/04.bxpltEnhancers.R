library(reshape2)
library(ggplot2)
library(ggsignif)
library(Biobase)
library(ggpubr)
### Loading complete eset
# Path to TCGA object
path <- "~/FANTOM6/TCGA/tcgaEset.rda"
load(path)
### Selecting genes of interest
DGEelnc <- c("CATG00000107122","CATG00000054627", "CATG00000039286", "ENSG00000231246", "ENSG00000255958")

### Subseting matrix and removing objs to clear RAM
exp <- exprs(TCGAeset)[DGEelnc,]
pheno <- pData(TCGAeset)
rm(TCGAeset)

### Removing cancers without normals
selTumors <- c("Bile Duct", "Bladder", "Breast", "Colorectal", "Esophagus", 
          "Head and Neck", "Kidney", "Liver", "Lung", "Prostate",
          "Stomach", "Thyroid", "Uterus")
keep <- pheno$gdc_cases.project.primary_site %in% selTumors
phenoSel <- pheno[keep,]
exp <- exp[,keep]
table(phenoSel$gdc_cases.project.primary_site)

### Removing metastatic and others tumor types
keep <- phenoSel$cgc_sample_sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
phenoSel <- phenoSel[keep,]
exp <- exp[,keep]
table(phenoSel$gdc_cases.project.primary_site)

### Creating df to plot
phenoSel <- cbind.data.frame(Site = phenoSel$gdc_cases.project.primary_site, Type = phenoSel$cgc_sample_sample_type)
df <- cbind.data.frame(t(exp), phenoSel)
df[,1:5] <- log2(df[,1:5]+2)
melt <- melt(df)

### Plotting fig
png(filename = "./figs/consensusBxplt_elnc.png", height = 3000, width = 6000, 
    res = 300)
plot1 <- ggplot(melt, aes(x = Site , y = value, fill = Type) ) +
        geom_boxplot(outlier.shape = NA) + ## geom_boxplot(outlier.shape = NA) +
        labs(title="Expression levels for consensus enhancers",
             x="", y="log2 counts") +
        stat_summary(fun.y=mean, geom="point", shape=124, size=5, 
                     color="White", position=position_dodge(width=0.75),) +
        facet_wrap(~ variable, scales = "fixed", ncol=5) +
        stat_compare_means(data = melt, method.args = 
                          list(formula = value ~ Type, 
                          group.by = c("Site", "variable"),
                          p.adjust.method = "BH"),
                          method = "t.test", label = "p.signif", 
                          symnum.args = list(cutpoints = 
                          c(0, 0.0001, 0.001, 0.01, 0.1, 1), symbols = c("****", 
                          "***", "**", "*", "ns"))) +
        coord_flip() +
        theme_light() +
        theme(plot.title = element_text(size = 25, face="bold", hjust = 0.5),
              axis.title = element_text(size=20),
              axis.text.x = element_text(size=18),
              axis.text.y = element_text(size=18),
              legend.key.size = unit(3, "lines"),
              legend.title = element_text(size=20, face="plain"),
              legend.text = element_text(size=15),
              strip.text = element_text(size=20)) 
plot1        
dev.off()
