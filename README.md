---

# Description
This repository hosts all the code used in the FC-R2 paper by Imada & Sanchez et al. 2019.
The code is organized in three folders. 
- Pre-Processing folder contain the scripts used to create and process the resource. 
- TCGA folder contain the scripts necessary for the analysis performed in the TCGA cohort, such as the DGE and Survival analysis (under subfolder survAnalysis). 
- GTEX folder contain the scripts necessary for all analysis performed in the GTEx cohort, this includes the specificity and expression levels analysis, biomarkers expression profiles across tissues and other validation steps. 

All data necessary to reproduce these analysis can be obtained in https://jhubiostatistics.shinyapps.io/recount/
## How to use

Data is available as a RangeSummarizedObject (RSE) which can be loaded in R enviroment using the SummarizedExperiment package.
Raw data is available as overall base coverage at gene level.

Data can be used as is, or scaled by the number of mapped reads and read length or AUC to get read-counts like data. For more information about processing data and generating recount compatible tracks see:

*Collado-Torres, Leonardo, Abhinav Nellore, and Andrew E. Jaffe. "recount workflow: Accessing over 70,000 human RNA-seq samples with Bioconductor." F1000Research 6 (2017).*
*Nellore, Abhinav, Leonardo Collado-Torres, Andrew E. Jaffe, José Alquicira-Hernández, Christopher Wilks, Jacob Pritt, James Morton, Jeffrey T. Leek, and Ben Langmead. "Rail-RNA: scalable analysis of RNA-seq splicing and coverage." Bioinformatics 33, no. 24 (2016): 4033-4040.*

```r
### Load required libraries
library(SummarizedExperiment)
library(recount)

### Download data
download_study('TCGA', type = 'rse-fc')

### Scale data
TCGA <- scale_counts(rse_fc, by="auc")
```

For huge datasets (e.g. TCGA or GTEX) you might want to subset only data relevant to your work (e.g. cancer type) to reduce the amount of computational resource used.

```r
### Load required libraries
library(SummarizedExperiment)
library(recount)

### Extract counts matrix
rseCounts <- assays(TCGA)$counts

### Extract feature information
grList <- rowRanges(GTEX)
fdata <- mcols(grList)


### Creating list of pheno data splited by cancer type
pheno <- as.data.frame(colData(TCGA))
phenoList <- split(pheno,pheno$gdc_cases.project.project_id)

### Check number of samples and projects available
table(pheno$gdc_cases.project.project_id)

### Create expression sets for a cancer type
x <- "TCGA_BRCA"
exp <- matrix[,rownames(phenoList[[x]])]
pheno <- phenoList[[x]]
### Create expression set for limma/voom
eset <- ExpressionSet(assayData <- exp,
                      phenoData = AnnotatedDataFrame(pheno),
                      featureData= AnnotatedDataFrame(fdata))
```

You can browse a study of interest at recount2 website or in R.
```
### Get table of all studies containing keywords
project_info <- abstract_search('HeLa cells')
### Select one of projects ID
download_study('SRP057804', type = 'rse-fc')
```
### Search for studies
## Citation
If you used any of the code or dataset available please cite:
*Imada, Eddie Luidy, et al. "Recounting the FANTOM Cage Associated Transcriptome." BioRxiv (2019): 659490.*
