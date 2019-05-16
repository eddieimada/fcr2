### Disjoining exons ranges by strand

### Setting WD & loading libraries
setwd("~/Dropbox/FANTOM6/GRangesAnnotation/")
library(dplyr)

### Load exon ranges and disjoin
load("./objs/exon_ranges.rda")
djn <- disjoin(exon_ranges)

### Remapping disjoined ranges to orignal exon ranges
remapping <- findOverlaps(exon_ranges, djn)

### Checking disjoined ranges mapped to more than one exon and subsetting 
### ambiguous and unique disjoined ranges
dup <- duplicated(subjectHits(remapping))
ambValuesIndex <- unique(subjectHits(remapping[dup]))
ambiguous <- remapping[subjectHits(remapping) %in% ambValuesIndex]
uniq <- remapping[!subjectHits(remapping) %in% ambValuesIndex]

### Converting GRanegs to df
uniq <- as.data.frame(uniq)
exon_ranges <- as.data.frame(exon_ranges)
ambiguous <- as.data.frame(ambiguous)

### Creating a DF for unique ranges annotated
uniqID <- exon_ranges[uniq$queryHits,]$id
uniq <- as.data.frame(djn[uniq$subjectHits])
uniqAnno <- cbind.data.frame(uniq,uniqID)

### Creating a DF for ambiguous ranges
ambiguous <- ambiguous[order(ambiguous$subjectHits),]
ambiguous <- cbind.data.frame(ambiguous, "id"=exon_ranges$id[ambiguous$queryHits])
ambiguous$id <- as.character(ambiguous$id)

### Create a list of lists, where each list contain all exons overlapping each
### disjoined range
subsets<-split(ambiguous, as.factor(ambiguous$subjectHits), drop=TRUE)


### Collapse each list into a single string joined by ";"
step <- 0
range_n <- length(subsets)
collapsedNamesList <- sapply(subsets, function(x){
    step <<- step + 1
    cat("\r", paste0("Processing range ", step, " of ", range_n))
    str <- paste(x$id, sep="", collapse = ";")
})


### Creating a DF of unique ambiguous ranges
djnAmbiguous <- as.data.frame(djn[unique(ambiguous$subjectHits)])

### Bind collapsed list which is the ID for each disjoined range
ambiguousAnno <- cbind.data.frame(djnAmbiguous, collapsedNamesList)

### Renaming columns
names(ambiguousAnno)[6] <- "ID"
names(uniqAnno)[6] <- "ID"

### Creating full dataset with unique and ambiguous ranges
allAnno <-  rbind.data.frame(uniqAnno,ambiguousAnno)
allAnno <- allAnno[order(allAnno$strand, allAnno$seqnames, allAnno$start),]
rownames(allAnno) <- 1:nrow(allAnno)

### Create GR object and save
djnStrandGR <- makeGRangesFromDataFrame(allAnno, keep.extra.columns = TRUE)
save(djnStrandGR, file="objs/disjoinedByStrand.rda")
