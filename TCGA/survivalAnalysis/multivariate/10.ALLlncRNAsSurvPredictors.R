###########################################################################
### Survival analysis by subtypes using intergenic-plncRNA
### 
### Diego F Sanchez
###
###########################################################################
### First steps
### Set working directory
setwd(".")
#/Users/diegosanchez/DropboxMech/FANTOM6/TCGA/survAnalysis/code
### Clean enviroments
rm(list=ls())

### Load libraries
library(Biobase)
require(ggplot2)
require(reshape2)

###########################################################################
### Summarizing results from Cox PH regression for lncRNAs types
###########################################################################
enhancerDF <- read.csv("../text/significantEnhancersByTCGAtypes.csv")
divergentPromDF <- read.csv("../text/divergentPlncRNA/significantDivergentPromByTCGAtypes.csv")
intergenicPromDF <- read.csv("../text/intergenicPlncRNA/significantIntergenicPromByTCGAtypes.csv")
enhancerDF <- enhancerDF[ order(enhancerDF$Types), ]
divergentPromDF <- divergentPromDF[ order(divergentPromDF$Types), ]
intergenicPromDF <- intergenicPromDF[ order(intergenicPromDF$Types), ]
lncRNAsurvPredictors <- cbind.data.frame(Types = as.character(enhancerDF$Types), 
                              Enhancers = enhancerDF$FDR....0.05,
                              'Divergent promoters' = divergentPromDF$FDR....0.05,
                              'Intergenic promoters' = intergenicPromDF$FDR....0.05,
                              Cases = enhancerDF$Cases,
                              Events = enhancerDF$Events,
                              Median.time = enhancerDF$Median.time)
lncRNAsurvPredictors <- lncRNAsurvPredictors[order(lncRNAsurvPredictors$Enhancers, decreasing = T), ]
lncRNAsurvPredictors$Types <-factor(lncRNAsurvPredictors$Types, 
                                    levels = lncRNAsurvPredictors$Types)

# Saving into significantlncRNAsByTCGAtypes.csv
write.csv(lncRNAsurvPredictors, file = "../text/significantlncRNAsByTCGAtypes.csv")


####
dfmelted <- melt(lncRNAsurvPredictors[,c('Types','Enhancers','Divergent promoters', 
                                         'Intergenic promoters')], id.vars = 1)

### Making plot for number of genes vs TCGA cancer types
jpeg(filename = "../figs/NumberOflncRNAsByTCGAtypes.jpg",  width = 4000, height = 3000, res = 300)
ggplot(data = dfmelted, aes(x = Types, y = value, group = variable)) +
  theme_classic() +
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
  geom_text(aes(label = value), hjust=-0.1, 
            angle = 90,
            size = 5,
            position = position_dodge(0.9)) +
  xlab("") + ylab("Number of lncRNAs") +
  theme(axis.text.x = element_text(color = "black", 
                                 size= 20, 
                                 angle=90, 
                                 vjust=.8, 
                                 hjust=0.8),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 20),
        legend.position= c(0.8,0.6),
        legend.title = element_blank()) +
  expand_limits(y = 4100) +
  ggtitle("")
dev.off()

