#################################################################
### Categorization of lncRNA's

#################################################################
### Clean
rm(list=ls())
##################################################################
### Load packages & data

require(Biobase)
# Path to GTEx object
path <- "~/FANTOM6/GTEX/gtexEset.rda"
load(path)

#################################################################
### Look at info of DHS and lncRNA

table(featureData(gtexEset)$CAT_DHS_type,featureData(gtexEset)$CAT_geneClass)

#################################################################
### Categorize lncRNA

##########
## Create intergenic p-lncRNA group

intergenicP<- gtexEset[featureData(gtexEset)$CAT_geneClass == "lncRNA_intergenic" &
                       featureData(gtexEset)$CAT_DHS_type == "DHS_promoter", ]

table(featureData(intergenicP)$CAT_DHS_type, featureData(intergenicP)$CAT_geneClass)

##########
## Create divergent p-lncRNA group 

divergentP <- gtexEset[featureData(gtexEset)$CAT_geneClass == "lncRNA_divergent" &
                         featureData(gtexEset)$CAT_DHS_type == "DHS_promoter", ]

table(featureData(divergentP)$CAT_DHS_type, featureData(divergentP)$CAT_geneClass)

##########
## Create e-lncRNA group

elncRNA <- gtexEset[featureData(gtexEset)$CAT_geneClass == "lncRNA_divergent" & featureData(gtexEset)$CAT_DHS_type == "DHS_enhancer" |
                  featureData(gtexEset)$CAT_geneClass == "lncRNA_intergenic" & featureData(gtexEset)$CAT_DHS_type == "DHS_enhancer" |
                  featureData(gtexEset)$CAT_geneClass == "lncRNA_antisense" & featureData(gtexEset)$CAT_DHS_type == "DHS_enhancer" |
                  featureData(gtexEset)$CAT_geneClass == "lncRNA_sense_intronic" & featureData(gtexEset)$CAT_DHS_type == "DHS_enhancer", ]

table(featureData(elncRNA)$CAT_DHS_type, featureData(elncRNA)$CAT_geneClass)

##########
## Create coding mRNA group

codingmRNA <- gtexEset[featureData(gtexEset)$CAT_geneClass == "coding_mRNA" & featureData(gtexEset)$CAT_DHS_type == "DHS_dyadic" |
                         featureData(gtexEset)$CAT_geneClass == "coding_mRNA" & featureData(gtexEset)$CAT_DHS_type == "DHS_enhancer" |
                         featureData(gtexEset)$CAT_geneClass == "coding_mRNA" & featureData(gtexEset)$CAT_DHS_type == "DHS_promoter" |
                         featureData(gtexEset)$CAT_geneClass == "coding_mRNA" & featureData(gtexEset)$CAT_DHS_type == "not_DHS", ]

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
