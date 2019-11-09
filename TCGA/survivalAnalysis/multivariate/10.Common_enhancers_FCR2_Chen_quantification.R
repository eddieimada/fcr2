###########################################################################
### FC-R2
### Genes in common Chen's paper and FC-R2
### Diego F Sanchez
###
###########################################################################


rm(list = ls())

### Mapping between Chen's and FANTOM6
Mapping_table <- read.csv("~/DropboxMech/FANTOM6/Manuscripts/FCR2/Pan_cancer/text/mapping_table.csv", 
                          header = T, stringsAsFactors = FALSE)
length(Mapping_table$fantom_cat)
### nnumber of genes in Chen's mapped to FCR2 genes (different types)
length(unique(Mapping_table$fantom_cat))

### Chen predictive enhancers coordinates
Chen_enhancers <- read.csv("~/DropboxMech/FANTOM6/Manuscripts/FCR2/Pan_cancer/text/ChenPredictiveEnhancers.csv",
                           header = T, stringsAsFactors = FALSE)
# number of prognostic genes in Chen's paper 
length(unique(Chen_enhancers$Prognostic_enhancer))

### FCR2 predictive enhancers
FCR2_enhancers <- read.csv("~/DropboxMech/FANTOM6/Manuscripts/FCR2/Pan_cancer/text/union_of_predictive_enhancers_through_tumors_list.csv",
                           header = T, stringsAsFactors = FALSE)
FCR2_enhancers <- FCR2_enhancers$x
# FCR2 number of predictive genes accross 13 cancer types.
length(FCR2_enhancers)
length(unique(FCR2_enhancers))

### Ehnancers in common (Chen's intersection FCR2)
Chen_to_FCR2 <- Mapping_table[Mapping_table$pan_cancer %in% unique(Chen_enhancers$Prognostic_enhancer), ]
# number of regions in common
length(unique(Chen_to_FCR2$pan_cancer))
# number of genes in common
length(unique(Chen_to_FCR2$fantom_cat))

     