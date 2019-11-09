###########################################################################
### Latex tables for predictive genes by RNA subtypes
### 
### Diego F Sanchez
###
###########################################################################
### First steps

### Set working directory
setwd(".")
#setwd(~/DropboxMech/FANTOM6/TCGA/survAnalysis/code)

### Clean enviroments
rm(list=ls())

### Load libraries
library(xtable)

###########################################################################
### load files
newRowNames <- c("Tumor type", "Non-significant", "FDR < 0.05", "Cases", "Events", "Median time")
mRNA_predictive <- read.csv(file="../text/mRNA/significantCodingmRNAByTCGAtypeslog.csv", 
                               stringsAsFactors = FALSE)
mRNA_predictive <- mRNA_predictive[,-1]
names(mRNA_predictive) <- newRowNames

elncRNA_predictive <- read.csv(file="../text/significantEnhancersByTCGAtypes_filtered_log.csv", 
                            stringsAsFactors = FALSE)
elncRNA_predictive <- elncRNA_predictive[,c(-1,-3)]
names(elncRNA_predictive) <- newRowNames

dplncRNA_predictive <- read.csv(file="../text/divergentPlncRNA/significantDivergentPromByTCGAtypes_log.csv", 
                               stringsAsFactors = FALSE)
dplncRNA_predictive <- dplncRNA_predictive[,-1]
names(dplncRNA_predictive) <- newRowNames

iplncRNA_predictive <- read.csv(file="../text/intergenicPlncRNA/significantIntergenicPromByTCGAtypes_log.csv", 
                                stringsAsFactors = FALSE)
iplncRNA_predictive <- iplncRNA_predictive[,-1]
names(iplncRNA_predictive) <- newRowNames


###########################################################################
### Creating latex files
### To print long tables prepare elements for longtable
addtorow <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(0)
addtorow$command  <- c(paste("\\hline \n", "\\endhead \n", "\\hline \n",
                             "{\\tiny Continued on next page} \n", "\\endfoot \n",
                             "\\endlastfoot \n", sep=""))
### Write mRNA latexTable
sink(file="../text/mRNA/significantmRNAByTCGAtypes_log.tex")
print(xtable(mRNA_predictive,
             caption = "Survival analysis using Cox proportional regression showing the number of mRNAs with prognostic value accross the 13 cancer types. \\textit{Non-significant} column indicates the number of genes with FDR greater than 0.05. \\textit{Cases} represents the number patients at the beginning of follow-up for each tumor type. \\textit{Events} is the number of death cases during follow up. \\textit{Median time} is given in days.",
             label = "tab:mRNA_predictive"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=FALSE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

### Write elncRNA latexTable
sink(file="../text/significantEnhancersByTCGAtypes_log.tex")
print(xtable(elncRNA_predictive,
             caption = "Survival analysis using Cox proportional regression showing the number of e-lncRNAs with prognostic value accross the 13 cancer types. \\textit{Non-significant} column indicates the number of genes with FDR greater than 0.05. \\textit{Cases} represents the number patients at the beginning of follow-up for each tumor type. \\textit{Events} is the number of death cases during follow up. \\textit{Median time} is given in days.",
             label = "tab:e-lncRNA_predictive"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=FALSE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

### Write dplncRNA latexTable
sink(file="../text/divergentPlncRNA/significantDivergentPromByTCGAtypes_log.tex")
print(xtable(dplncRNA_predictive,
             caption = "Survival analysis using Cox proportional regression showing the number of dp-lncRNAs with prognostic value accross the 13 cancer types. \\textit{Non-significant} column indicates the number of genes with FDR greater than 0.05. \\textit{Cases} represents the number patients at the beginning of follow-up for each tumor type. \\textit{Events} is the number of death cases during follow up. \\textit{Median time} is given in days.",
             label = "tab:dp-lncRNA_predictive"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=FALSE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

### Write dplncRNA latexTable
sink(file="../text/intergenicPlncRNA/significantIntergenicPromByTCGAtypes_log.tex")
print(xtable(iplncRNA_predictive,
             caption = "Survival analysis using Cox proportional regression showing the number of ip-lncRNAs with prognostic value accross the 13 cancer types. \\textit{Non-significant} column indicates the number of genes with FDR greater than 0.05. \\textit{Cases} represents the number patients at the beginning of follow-up for each tumor type. \\textit{Events} is the number of death cases during follow up. \\textit{Median time} is given in days.",
             label = "tab:ip-lncRNA_predictive"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=FALSE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

