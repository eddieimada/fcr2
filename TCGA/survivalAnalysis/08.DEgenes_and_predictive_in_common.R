###########################################################################
### Exploring predictive genes through RNA subtypes
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
#library(Biobase)
library(xtable)

### Xtable Latex preamble
# To print long tables prepare elements for longtable
addtorow <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(0)
addtorow$command  <- c(paste("\\hline \n", "\\endhead \n", "\\hline \n",
                             "{\\tiny Continued on next page} \n", "\\endfoot \n",
                             "\\endlastfoot \n", sep=""))

### Load files
# DGE files
mRNA_DGE_dataframe <- read.csv(file="../../textFinal/mRNA_geneIntersection.csv", 
                                  stringsAsFactors = FALSE)

elncRNA_DGE_dataframe <- read.csv(file="../../textFinal/elncRNA_geneIntersection.csv", 
                               stringsAsFactors = FALSE)

dplncRNA_DGE_dataframe <- read.csv(file="../../textFinal/dProm_geneIntersection.csv", 
                               stringsAsFactors = FALSE)

iplncRNA_DGE_dataframe <- read.csv(file="../../textFinal/iProm_geneIntersection.csv", 
                               stringsAsFactors = FALSE)
# Predictive files
mRNA_predictive <- read.csv(file="../text/mRNA/union_of_predictive_CodingmRNA_through_tumors_list.csv", 
                               stringsAsFactors = FALSE)
mRNA_predictive <- mRNA_predictive$x

elncRNA_predictive <- read.csv(file="../text/union_of_predictive_enhancers_through_tumors_list.csv", 
                            stringsAsFactors = FALSE)
elncRNA_predictive <- elncRNA_predictive$x

dplncRNA_predictive <- read.csv(file="../text/divergentPlncRNA/union_of_predictive_dplncRNAs_through_tumors_list.csv", 
                               stringsAsFactors = FALSE)
dplncRNA_predictive <- dplncRNA_predictive$x

iplncRNA_predictive <- read.csv(file="../text/intergenicPlncRNA/union_of_predictive_iplncRNA_through_tumors_list.csv", 
                               stringsAsFactors = FALSE)
iplncRNA_predictive <- iplncRNA_predictive$x

###########################################################################
### Looking for DE genes that are prognostic
mRNA_names <- mRNA_DGE_dataframe$geneID[mRNA_DGE_dataframe$geneID %in% mRNA_predictive]
length(mRNA_names)

elncRNA_names <- elncRNA_DGE_dataframe$geneID[elncRNA_DGE_dataframe$geneID %in% elncRNA_predictive]
length(elncRNA_names)

dplncRNA_names <- dplncRNA_DGE_dataframe$geneID[dplncRNA_DGE_dataframe$geneID %in% dplncRNA_predictive]
length(dplncRNA_names)

iplncRNA_names <- iplncRNA_DGE_dataframe$geneID[iplncRNA_DGE_dataframe$geneID %in% iplncRNA_predictive]
length(iplncRNA_names)

###########################################################################
### mRNA predictive accross tumors
mRNA_files <- list.files("~/DropboxMech/FANTOM6/TCGA/survAnalysis/text/mRNA", recursive = TRUE,
                    pattern="Sig.+csv$", full.names = TRUE)
### Read files
mRNA <- lapply(mRNA_files, function(x) {
  read.csv(x, stringsAsFactors = FALSE)}) 
#mRNA <- unlist(mRNA, recursive = FALSE)

### Add names
nmsT <- gsub("\\.csv", "", gsub(".+mRNA_", "", mRNA_files))

names(mRNA) <- nmsT

# Removing uninformative DF (bile and esophagus)
mRNA <- mRNA[c(-1, -5)]

# # Extracting HR by types
# predictive <- lapply(mRNA, function(x){
#   x[x$X %in% mRNA_names, ]
# })

#Extracting HR by types
predictive <- lapply(mRNA, function(x){
  rownames(x) <- x$X
  HR <- x[mRNA_names, "HR" ]
  names(HR) <- mRNA_names
  HR
})
# bind vectors from list
mRNA_df <- do.call("cbind", predictive)
mRNA_df <- apply(mRNA_df,2, function(x){
  # rename
  ifelse(x <1, "Good", "Bad")
})
mRNA_df <- apply(mRNA_df,2, function(x){
  # rename
  ifelse(is.na(x), "N.P.", x)
})

### Write mRNA latexTable
sink(file="../text/mRNA/DGE_Prognostic_mRNAByTCGAtypes.tex")
print(xtable(mRNA_df,
             caption = "Differentially expressed mRNA genes with prognostic value accross cancer types. \\textit{Good} indicates the Cox HR $<$ 1. \\textit{Bad} represents Cox HR $>$ 1. \\textit{N.P.} refers that the given gene is non-prognostic on this cancer type",
             label = "tab:DGE_mRNA_predictive"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=TRUE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

### Write mRNA Table
write.csv(mRNA_df, file="../text/DGE_Prognostic_mRNAByTCGAtypes.csv")

###########################################################################
### e-lncRNA predictive accross tumors
elncRNA_files <- list.files("~/DropboxMech/FANTOM6/TCGA/survAnalysis/text/", recursive = FALSE,
                         pattern="Sig.+csv$", full.names = TRUE)
### Read files
elncRNA <- lapply(elncRNA_files, function(x) {
  read.csv(x, stringsAsFactors = FALSE)}) 
#mRNA <- unlist(mRNA, recursive = FALSE)

### Add names
nmsT <- gsub("\\.csv", "", gsub(".+RNA_", "", elncRNA_files))
names(elncRNA) <- nmsT

# Removing uninformative DF (bile and esophagus)
elncRNA <- elncRNA[c(-1, -5)]

#Extracting HR by types
predictive <- lapply(elncRNA, function(x){
  rownames(x) <- x$X
  HR <- x[elncRNA_names, "HR" ]
  names(HR) <- elncRNA_names
  HR
})
# bind vectors from list
elncRNA_df <- do.call("cbind", predictive)
elncRNA_df <- apply(elncRNA_df,2, function(x){
  # rename
  ifelse(x <1, "Good", "Bad")
})
elncRNA_df <- apply(elncRNA_df,2, function(x){
  # rename
  ifelse(is.na(x), "N.P.", x)
})

### Write elncRNA latexTable
sink(file="../text/DGE_Prognostic_elncRNAByTCGAtypes.tex")
print(xtable(elncRNA_df,
             caption = "Differentially expressed e-lncRNA genes with prognostic value accross cancer types. \\textit{Good} indicates the Cox HR $<$ 1. \\textit{Bad} represents Cox HR $>$ 1. \\textit{N.P.} refers that the given gene is non-prognostic on this cancer type",
             label = "tab:DGE_elncRNA_predictive"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=TRUE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

### Write elncRNA Table
write.csv(elncRNA_df, file="../text/DGE_Prognostic_elncRNAByTCGAtypes.csv")

###########################################################################
### dp-lncRNA predictive accross tumors
dplncRNA_files <- list.files("~/DropboxMech/FANTOM6/TCGA/survAnalysis/text/divergentPlncRNA/", recursive = FALSE,
                            pattern="Sig.+csv$", full.names = TRUE)
### Read files
dplncRNA <- lapply(dplncRNA_files, function(x) {
  read.csv(x, stringsAsFactors = FALSE)}) 

### Add names
nmsT <- gsub("\\.csv", "", gsub(".+RNA_", "", dplncRNA_files))
names(dplncRNA) <- nmsT

# Removing uninformative DF (bile, esophagus, and Head and Neck)
dplncRNA <- dplncRNA[c(-1, -5, -6)]

#Extracting HR by types
predictive <- lapply(dplncRNA, function(x){
  rownames(x) <- x$X
  HR <- x[dplncRNA_names, "HR" ]
  names(HR) <- dplncRNA_names
  HR
})
# bind vectors from list
dplncRNA_df <- do.call("cbind", predictive)
dplncRNA_df <- apply(dplncRNA_df,2, function(x){
  # rename
  ifelse(x <1, "Good", "Bad")
})
dplncRNA_df <- apply(dplncRNA_df,2, function(x){
  # rename
  ifelse(is.na(x), "N.P.", x)
})

### Write dplncRNA latexTable
sink(file="../text/divergentPlncRNA/DGE_Prognostic_dplncRNAByTCGAtypes.tex")
print(xtable(dplncRNA_df,
             caption = "Differentially expressed dp-lncRNA genes with prognostic value accross cancer types. \\textit{Good} indicates the Cox HR $<$ 1. \\textit{Bad} represents Cox HR $>$ 1. \\textit{N.P.} refers that the given gene is non-prognostic on this cancer type",
             label = "tab:DGE_dplncRNA_predictive"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=TRUE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

### Write dplncRNA Table
write.csv(dplncRNA_df, file="../text/divergentPlncRNA/DGE_Prognostic_dplncRNAByTCGAtypes.csv")

###########################################################################
### ip-lncRNA predictive accross tumors
iplncRNA_files <- list.files("~/DropboxMech/FANTOM6/TCGA/survAnalysis/text/intergenicPlncRNA/", recursive = FALSE,
                             pattern="Sig.+csv$", full.names = TRUE)
### Read files
iplncRNA <- lapply(iplncRNA_files, function(x) {
  read.csv(x, stringsAsFactors = FALSE)}) 

### Add names
nmsT <- gsub("\\.csv", "", gsub(".+RNA_", "", iplncRNA_files))
names(iplncRNA) <- nmsT

# Removing uninformative DF (bile and esophagus)
iplncRNA <- iplncRNA[c(-1, -5)]

#Extracting HR by types
predictive <- lapply(iplncRNA, function(x){
  rownames(x) <- x$X
  HR <- x[iplncRNA_names, "HR" ]
  names(HR) <- iplncRNA_names
  HR
})
# bind vectors from list
iplncRNA_df <- do.call("cbind", predictive)
iplncRNA_df <- apply(iplncRNA_df,2, function(x){
  # rename
  ifelse(x <1, "Good", "Bad")
})
iplncRNA_df <- apply(iplncRNA_df,2, function(x){
  # rename
  ifelse(is.na(x), "N.P.", x)
})

### Write iplncRNA latexTable
sink(file="../text/intergenicPlncRNA/DGE_Prognostic_iplncRNAByTCGAtypes.tex")
print(xtable(iplncRNA_df,
             caption = "Differentially expressed ip-lncRNA genes with prognostic value accross cancer types. \\textit{Good} indicates the Cox HR $<$ 1. \\textit{Bad} represents Cox HR $>$ 1. \\textit{N.P.} refers that the given gene is non-prognostic on this cancer type",
             label = "tab:DGE_iplncRNA_predictive"), #align = rep("l", ncol(tmp)))
      type="latex", include.rownames=TRUE, table.placement="H",
      caption.placement="top", size="scriptsize", NA.string="N.A.",
      tabular.environment='longtable', floating=FALSE,
      add.to.row = addtorow,  hline.after=c(-1))
sink(file=NULL)

### Write iplncRNA Table
write.csv(iplncRNA_df, file="../text/intergenicPlncRNA/DGE_Prognostic_iplncRNAByTCGAtypes.csv")


###########################################################################
### Session information and clean quit

### Session
sessionInfo()

### Clean
rm(list=ls())

### Quit
q("no")
