####################################################################
### Code for subsetting by cancer type
####################################################################

### Clean and set wd
rm(list=ls())
library(Biobase)
### Setting WD and loading data
# Path to object
path <- "~FANTOM6/TCGA/codingmRNA.rda"
load(path)


### Creating list of pheno data splited by cancer type
phenoList <- split(pData(codingmRNA),pData(codingmRNA)$gdc_cases.project.primary_site)

### Removing NA's
NAindex <- which(is.na(rownames(exprs(codingmRNA))))
matrix <- exprs(codingmRNA)[-NAindex,]
fdata <- fData(codingmRNA)[-NAindex,]

### Storing list names
nm <- names(phenoList)

### Create separate eset's for each cancer type
list <- lapply(1:length(phenoList), function(x){
    exp <- matrix[,rownames(phenoList[[x]])]
    pheno <- phenoList[[x]]
    eset <- ExpressionSet(assayData <- exp,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(fdata))
    obj <- paste0(nm[[x]],"CodingmRNA")
    obj <- gsub(" ", "", obj)
    assign(obj,eset, envir = .GlobalEnv)
})


### Get all newly created objs and save them
objs <- ls()[grepl("CodingmRNA",ls())]
objs <- objs[!grepl("^CodingmRNA", objs)]
sapply(objs, function(x){
    save(list=x, file = paste0("./objs/", x, ".rda"))    
})    


############################################################################
### Clean and set wd
rm(list=ls())

### Loading data
# Path to object
path <- "~FANTOM6/TCGA/intergenicP.rda"
load(path)

### Creating list of pheno data splited by cancer type
phenoList <- split(pData(intergenicP),pData(intergenicP)$gdc_cases.project.primary_site)

### Removing NA's
NAindex <- which(is.na(rownames(exprs(intergenicP))))
matrix <- exprs(intergenicP)[-NAindex,]
fdata <- fData(intergenicP)[-NAindex,]

### Storing list names
nm <- names(phenoList)

### Create separate eset's for each cancer type
list <- lapply(1:length(phenoList), function(x){
    exp <- matrix[,rownames(phenoList[[x]])]
    pheno <- phenoList[[x]]
    eset <- ExpressionSet(assayData <- exp,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(fdata))
    obj <- paste0(nm[[x]],"intergenicP")
    obj <- gsub(" ", "", obj)
    assign(obj,eset, envir = .GlobalEnv)
})


### Get all newly created objs and save them
objs <- ls()[grepl("intergenicP",ls())]
objs <- objs[!grepl("^intergenicP", objs)]
sapply(objs, function(x){
    save(list=x, file = paste0("./objs/", x, ".rda"))    
})    


############################################################################
### Clean and set wd
rm(list=ls())

### Loading data
# Path to object
path <- "~FANTOM6/TCGA/elncRNA.rda"
load(path)


### Creating list of pheno data splited by cancer type
phenoList <- split(pData(elncRNA),pData(elncRNA)$gdc_cases.project.primary_site)

### Removing NA's
NAindex <- which(is.na(rownames(exprs(elncRNA))))
matrix <- exprs(elncRNA)[-NAindex,]
fdata <- fData(elncRNA)[-NAindex,]

### Storing list names
nm <- names(phenoList)

### Create separate eset's for each cancer type
list <- lapply(1:length(phenoList), function(x){
    exp <- matrix[,rownames(phenoList[[x]])]
    pheno <- phenoList[[x]]
    eset <- ExpressionSet(assayData <- exp,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(fdata))
    obj <- paste0(nm[[x]],"elncRNA")
    obj <- gsub(" ", "", obj)
    assign(obj,eset, envir = .GlobalEnv)
})


### Get all newly created objs and save them
objs <- ls()[grepl("elncRNA",ls())]
objs <- objs[!grepl("^elncRNA", objs)]
sapply(objs, function(x){
    save(list=x, file = paste0("./objs/", x, ".rda"))    
})    

############################################################################
### Clean and set wd
rm(list=ls())

### Loading data
# Path to object
path <- "~FANTOM6/TCGA/divergentP.rda"
load(path)


### Creating list of pheno data splited by cancer type
phenoList <- split(pData(divergentP),pData(divergentP)$gdc_cases.project.primary_site)

### Removing NA's
NAindex <- which(is.na(rownames(exprs(divergentP))))
matrix <- exprs(divergentP)[-NAindex,]
fdata <- fData(divergentP)[-NAindex,]

### Storing list names
nm <- names(phenoList)

### Create separate eset's for each cancer type
list <- lapply(1:length(phenoList), function(x){
    exp <- matrix[,rownames(phenoList[[x]])]
    pheno <- phenoList[[x]]
    eset <- ExpressionSet(assayData <- exp,
                          phenoData = AnnotatedDataFrame(pheno),
                          featureData= AnnotatedDataFrame(fdata))
    obj <- paste0(nm[[x]],"divergentP")
    obj <- gsub(" ", "", obj)
    assign(obj,eset, envir = .GlobalEnv)
})


### Get all newly created objs and save them
objs <- ls()[grepl("divergentP",ls())]
objs <- objs[!grepl("^divergentP", objs)]
sapply(objs, function(x){
    save(list=x, file = paste0("./objs/", x, ".rda"))    
})    

