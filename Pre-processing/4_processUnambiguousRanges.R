### Processing genes to ranges mapping, keeping only unambiguous ranges at 
### gene level

# Set wd and load objs
setwd("/home/elimada/Dropbox/FANTOM6/GRangesAnnotation/")
load("./objs/disjoinedByunstrand.rda")

rm(list=ls())
### Load mapping data
df <- read.csv("./data/full_table.tsv.gz", sep = "\t")

### Create ranges IDs without hash, and factors as integers (much faster processing)

df$id <- gsub("_.+", "", df$ranges)
df$ranges_num <- as.numeric(factor(df$id))
df$genes_num <- as.numeric(df$genes)

### Order df
df <- df[order(df$ranges_num),]

### Subset genes by ranges, and keep only those with only one ID
subsets <- split(df$genes_num, df$ranges_num, drop=TRUE)

keep <- lapply(subsets, function (x) {
        all(x == x[1])
})

### Expand keep to table size
keepExt <- sapply(df$ranges_num, function (x, y){
        y[[x]]
}, y = keep)

### Subset unambiguous ranges
unambiguous <- df[keepExt,]
### Remove ranges contributing to same gene more than once
uniqRanges <- unambiguous[!duplicated(unambiguous$ranges),]

### Getting width of each range
unambGR <- djnUnstrandGR[djnUnstrandGR$ID %in% uniqRanges$ranges]
rangeWidth <- width(unambGR)
names(rangeWidth) <- unambGR$ID
rangeWidth <- rangeWidth[as.vector(uniqRanges$ranges)]
uniqRanges <- cbind.data.frame(uniqRanges, rangeWidth)
uniqRanges <- droplevels(uniqRanges)

### Getting total number of bases covered for each gene
GeneCovFromRanges <- tapply(uniqRanges$rangeWidth, uniqRanges$genes, FUN=sum)
GeneCovFromRanges <- GeneCovFromRanges[uniqRanges$genes]
uniqRanges <- cbind.data.frame(uniqRanges, GeneCovFromRanges)


### Subset by strand
uniqRangesPos <- uniqRanges[uniqRanges$strand == "+",]
uniqRangesPos <-  droplevels(uniqRangesPos)
uniqRangesNeg <- uniqRanges[uniqRanges$strand == "-",]
uniqRangesNeg <-  droplevels(uniqRangesNeg)

### Save objs
save(uniqRangesNeg, file = "objs/uniqRangesNeg.rda")
save(uniqRangesPos, file = "objs/uniqRangesPos.rda")
