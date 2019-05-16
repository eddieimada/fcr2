### Script to create exons GRanges from CAT FANTOM data
library(GenomicRanges)
### setting wd
setwd("~/Dropbox/FANTOM6/GRangesAnnotation/")

### Reading transcripts annotation
transcripts <-
    read.delim("data/F6_CAT.transcript.bed.gz",
               sep = "\t",
               header = FALSE)

### Naming columns
names(transcripts) <-
    c(
        "chr",
        "start",
        "stop",
        "trnscID",
        "NA",
        "strand",
        "start",
        "stop",
        "NA2",
        "exon_num",
        "width",
        "baseFromStart"
    )

### Removing transcripts with the exact coordinates
#head(transcripts[order(rev(transcripts$chr),rev(transcripts$start), transcripts$exon_num, decreasing = TRUE), ])
#transcripts <- transcripts[!duplicated(paste0(transcripts$chr,transcripts$start,transcripts$stop)),]

### Changing vectors class
transcripts$start <- as.numeric(transcripts$start)
transcripts$stop <- as.numeric(transcripts$stop)

###  Variable to keep track of processing
step <- 0
transcripts_n <- nrow(transcripts)
### For each transcript, get exons information and save coordinates to a file
apply(transcripts, 1, function(x) {
    ### Getting exons positions and width for each transcript
    exon_starts <- as.numeric(unlist(strsplit(x[[12]], ",")))
    widths <- as.numeric(unlist(strsplit(x[[11]], ",")))
    ### for each transcript create start, stop, ID and strand info
    for (i in 1:x[[10]]) {
        exon_start <- as.numeric(x[[2]]) + exon_starts[i]
        exon_end <- exon_start + widths[[i]]
        ID <-  paste0(x[[4]], ".E", i)
        p <-
            paste(x[[1]], exon_start, exon_end, ID, paste0(x[[6]], "\n"), sep = "\t")
        cat(p, file = "data/exon_ranges.tsv", append = TRUE)
    }
    ### Keeping track of progress
    step <<- step + 1
    cat("\r",
        paste0("Processing transcript ", step, " of ", transcripts_n))
})

### Reading exon coordinates from file (it's faster to save to a file and read than assigning directly to object)
exon_ranges <- read.table("data/exon_ranges.tsv.gz", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
names(exon_ranges) <-  c("chr", "start", "stop", "id", "strand")

### Reading gene coordinates from file
gene_ranges <- read.table("data/F6_CAT.gene.bed.gz", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
names(gene_ranges) <-  c("chr", "start", "stop", "id", "n_transcripts", "strand")

### Creating GRanges objects from DF
transcripts_ranges <- makeGRangesFromDataFrame(transcripts[,-c(5,7,8,9)], keep.extra.columns = TRUE)
save(transcripts_ranges, file = "objs/transcripts_ranges.rda")
exon_ranges <- makeGRangesFromDataFrame(exon_ranges, keep.extra.columns = TRUE)
save(exon_ranges, file = "objs/exon_ranges.rda")
gene_ranges <- makeGRangesFromDataFrame(gene_ranges, keep.extra.columns = TRUE)
save(gene_ranges, file = "objs/gene_ranges.rda")



exon_rangesNR <-exon_ranges[!duplicated(paste0(seqnames(exon_ranges), start(exon_ranges), end(exon_ranges), strand(exon_ranges))),]
save(exon_rangesNR, file = "objs/exon_rangesNR.rda")

gene_rangesNR <- gene_ranges[!duplicated(paste0(seqnames(gene_ranges), start(gene_ranges), end(gene_ranges))),]
save(gene_rangesNR, file = "objs/gene_rangesNR.rda")

### Finding overlaps between "genes" ranges and exons ranges
hits <- findOverlaps(gene_rangesNR, exon_rangesNR)

### Keeping only mappings where exon in entirely contained inside gene range
overlaps <- pintersect(gene_ranges[queryHits(hits)], exon_ranges[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(exon_ranges[subjectHits(hits)])
hits <- hits[percentOverlap == 1]

### Creating DF with gene to exons mappings
geneToExon <- cbind.data.frame("genes"=gene_ranges$id[queryHits(hits)], "exons"=exon_ranges$id[subjectHits(hits)])
save(geneToExon, file = "objs/geneToExonMapping.rda")

### Finding overlaps between exons
hitsE2E <- findOverlaps(exon_rangesNR, exon_rangesNR)

### Removing self hits
hitsE2E <- hitsE2E[!queryHits(hitsE2E) == subjectHits(hitsE2E)]
save(hitsE2E, file = "objs/E2Eoverlaps.rda")
### Getting non-overlapping exons
nonOverlapExons <- exon_rangesNR[-unique(queryHits(hitsE2E))]
save(nonOverlapExons, file = "objs/nonOverlappingExons.rda")








