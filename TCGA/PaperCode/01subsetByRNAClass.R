####################################################################
### Code for subsetting original matrix into RNA class type matrices
####################################################################
### Clean and set wd
rm(list=ls())


##################################################################
### Load packages & data

require(Biobase)
# Path to TCGA object
path <- "~/FANTOM6/TCGA/tcgaEset.rda"
load(path)

#################################################################
### Look at info of DHS and lncRNA

table(featureData(TCGAeset)$CAT_DHS_type,featureData(TCGAeset)$CAT_geneClass)

#################################################################
### Categorize lncRNA

##########
## Create intergenic p-lncRNA group

intergenicP<- TCGAeset[featureData(TCGAeset)$CAT_geneClass == "lncRNA_intergenic" &
                           featureData(TCGAeset)$CAT_DHS_type == "DHS_promoter", ]

table(featureData(intergenicP)$CAT_DHS_type, featureData(intergenicP)$CAT_geneClass)

##########
## Create divergent p-lncRNA group 

divergentP <- TCGAeset[featureData(TCGAeset)$CAT_geneClass == "lncRNA_divergent" &
                           featureData(TCGAeset)$CAT_DHS_type == "DHS_promoter", ]

table(featureData(divergentP)$CAT_DHS_type, featureData(divergentP)$CAT_geneClass)

##########
## Create e-lncRNA group

elncRNA <- TCGAeset[featureData(TCGAeset)$CAT_geneClass == "lncRNA_divergent" & featureData(TCGAeset)$CAT_DHS_type == "DHS_enhancer" |
                        featureData(TCGAeset)$CAT_geneClass == "lncRNA_intergenic" & featureData(TCGAeset)$CAT_DHS_type == "DHS_enhancer" |
                        featureData(TCGAeset)$CAT_geneClass == "lncRNA_antisense" & featureData(TCGAeset)$CAT_DHS_type == "DHS_enhancer" |
                        featureData(TCGAeset)$CAT_geneClass == "lncRNA_sense_intronic" & featureData(TCGAeset)$CAT_DHS_type == "DHS_enhancer", ]

table(featureData(elncRNA)$CAT_DHS_type, featureData(elncRNA)$CAT_geneClass)

##########
## Create coding mRNA group

codingmRNA <- TCGAeset[featureData(TCGAeset)$CAT_geneClass == "coding_mRNA" & featureData(TCGAeset)$CAT_DHS_type == "DHS_dyadic" |
                           featureData(TCGAeset)$CAT_geneClass == "coding_mRNA" & featureData(TCGAeset)$CAT_DHS_type == "DHS_enhancer" |
                           featureData(TCGAeset)$CAT_geneClass == "coding_mRNA" & featureData(TCGAeset)$CAT_DHS_type == "DHS_promoter" |
                           featureData(TCGAeset)$CAT_geneClass == "coding_mRNA" & featureData(TCGAeset)$CAT_DHS_type == "not_DHS", ]

table(featureData(codingmRNA)$CAT_DHS_type, featureData(codingmRNA)$CAT_geneClass)


#################################################################
## Save objects 

save(intergenicP, file = "objs/intergenicP.rda")
save(divergentP, file = "objs/divergentP.rda")
save(elncRNA, file = "objs/elncRNA.rda")
save(codingmRNA, file = "objs/codingmRNA.rda") 

#################################################################
### Session information and clean quit

### Session
sessionInfo()

### Clean
rm(list=ls())

### Quit
q("no")
