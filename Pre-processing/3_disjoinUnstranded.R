### Disjoining exons ranges by both strands

### Setting WD & loading libraries
setwd("~/Dropbox/FANTOM6/GRangesAnnotation/")
library(dplyr)

### Load exon ranges and disjoin
load("./objs/disjoinedByStrand.rda")
djn <- disjoin(djnStrandGR, ignore.strand = TRUE)

### Remapping disjoined ranges to orignal exon ranges
remapping <- findOverlaps(djnStrandGR, djn, ignore.strand=TRUE)

### Checking disjoined ranges mapped to more than one exon and subsetting 
### ambiguous and unique disjoined ranges
dup <- duplicated(subjectHits(remapping))
ambValuesIndex <- unique(subjectHits(remapping[dup]))
ambiguous <- remapping[subjectHits(remapping) %in% ambValuesIndex]
uniq <- remapping[!subjectHits(remapping) %in% ambValuesIndex]

### Converting GRanegs to df
uniq <- as.data.frame(uniq)
djnStrandGR <- as.data.frame(djnStrandGR)
uniqID <- djnStrandGR[uniq$queryHits,]$ID
uniqStrand <- djnStrandGR[uniq$queryHits,]$strand

### Creating a DF for unique ranges annotated
uniq <- as.data.frame(djn[uniq$subjectHits])
uniqAnno <- cbind.data.frame(uniq,uniqID)
uniqAnno$strand <- uniqStrand

### Creating a DF for ambiguous ranges
ambiguous <- as.data.frame(ambiguous)
ambiguous <- ambiguous[order(ambiguous$subjectHits),]
ambiguous <- cbind.data.frame(ambiguous, "ID"=djnStrandGR$ID[ambiguous$queryHits])
ambiguous$ID <- as.character(ambiguous$ID)

### Create a list of lists, where each list contain all exons overlapping each
### disjoined range
subsets <- split(ambiguous, as.factor(ambiguous$subjectHits), drop=TRUE)

### Collapse each list into a single string joined by ";"
step <- 0
range_n <- length(subsets)
collapsedNamesList <- sapply(subsets, function(x){
    step <<- step + 1
    cat("\r", paste0("Processing range ", step, " of ", range_n))
    str <- paste(x$ID, sep="", collapse = ";")
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
djnUnstrandGR <- makeGRangesFromDataFrame(allAnno, keep.extra.columns = TRUE)
save(djnUnstrandGR, file="objs/disjoinedByunstrand.rda")
