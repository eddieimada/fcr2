### Plot intersection boxplots
##############################

### Load Libraries
library(ggplot2)
library(reshape2)
library(Biobase)
library(gtable)
library(grid)
# setwd
setwd("~/Dropbox (MechPred)/FANTOM6/TCGA/")
rm(list=ls())
# load objs
load("objs/tcgaEset.rda")
# load results from survival analysis
load("objs/BHintersectList.rda")
# extract data
exp <- exprs(TCGAeset)
pheno <- pData(TCGAeset)
feat <- fData(TCGAeset)
rm(TCGAeset)
gc()
mapName <- setNames(feat$geneID,feat$geneName)
# select only tumors with > 10 normals
site <- pheno$gdc_cases.project.primary_site
sitekeep <- site %in% c("Bladder", "Bile", "Breast", "Colorectal", "Esophagus", "Head and Neck", "Kidney", "Liver", "Lung", "Prostate", "Stomach", "Thyroid", "Uterus")
# select only normal and tumors
type <- pheno$gdc_cases.samples.sample_type
typekeep <- type %in% c("Primary Tumor", "Solid Tissue Normal")
# select only normal and tumors on cancers with > 10 normals
keep <- sitekeep & typekeep
# subset
exp <- exp[,keep]
pheno <- pheno[keep,]
# Get tumor & cancer type
site <- pheno$gdc_cases.project.primary_site
type <- pheno$gdc_cases.samples.sample_type

up <- sapply(BHlistUP, function(x) x[1])
dn <- sapply(BHlistDN, function(x) x[1])
up <- up[c(4,1,3,2)]
dn <- dn[c(4,1,3,2)]

sel <- c(up,dn)
sel <- c(rbind(up,dn))
ids <- mapName[sel]

names(ids)[c(6,8)] <- c("MIR99AHG","GABARAPL1-AS1")
names(ids)[3:5] <- c("SLBP-DT","CYP2U1-AS1", "LINC01331")
expSel <- exp[ids,]
expSel <- t(expSel)
colnames(expSel) <- names(ids)
df <- cbind.data.frame(expSel+1,site, type)
df <- melt(df)
df$class <- NA
df$class[df$variable %in% ids[1:2]] <- "mRNA"
df$class[df$variable %in% ids[3:4]] <- "d-lncRNA"
df$class[df$variable %in% ids[5:6]] <- "i-lncRNA"
df$class[df$variable %in% ids[7:8]] <- "e-lncRNA"

p1 <- ggplot(df, aes(x=site, y=log2(value), fill = type)) +
    geom_boxplot(outlier.size = 0.25) +
    theme_bw() +
    facet_wrap(~variable, nrow = 1, scales = "free_x") +
    theme(legend.position="bottom", 
          legend.title = element_blank(), 
          text = element_text(size = 14),
          strip.text = element_text(colour = "white", face = "italic"),
          axis.text.x=element_text(size=14, angle=60, hjust=1),
          axis.text.y=element_text(size=18),
          plot.title = element_text(hjust=0.5)) +
    xlab("") +
    #ggtitle(i) +
    coord_flip() +
    ylab("log2 expression") +
    scale_fill_manual(values = c("deeppink1", "deepskyblue1"))
#Extract text grobs
g <- ggplot_gtable(ggplot_build(p1))
# grab strips
strip_t <- which(grepl('strip-t', g$layout$name))
# set colors
fills <- c("red", "red","purple","purple", "cyan4","cyan4", "darkgreen", "darkgreen")
### Replace grobs
k <- 1
for (i in strip_t) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
}
jpeg(filename = "figs/intersection.jpg", width = 6000, height = 3000, res =300)
grid.draw(g)
dev.off()

fills <- c("firebrick4", "firebrick1","darkorchid4","darkorchid3", "cyan4","cyan3", "green4", "green3")
### Replace grobs
k <- 1
for (i in strip_t) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
}
png(filename = "figs/intersection-paired.png", width = 6000, height = 3000, res =300)
grid.draw(g)
dev.off()
