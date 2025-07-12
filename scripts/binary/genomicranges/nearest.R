
source("lib/helpers.R")
source("scripts/reading/genomicranges.R")

library(GenomicRanges)
library(GenomeInfoDb)

files <- get_binary_arguments()
## ------------------------------------------------------------------
## 1. Read & harmonise
## ------------------------------------------------------------------
annotations <- read_genomic_file(files$annotation_file)
reads       <- read_genomic_file(files$read_file)

## a) keep only chromosomes used by *both* objects
common <- intersect(seqlevelsInUse(reads), seqlevelsInUse(annotations))

reads       <- keepSeqlevels(reads,       common, pruning.mode = "tidy")
annotations <- keepSeqlevels(annotations, common, pruning.mode = "tidy")

hits <- distanceToNearest(reads, annotations, ignore.strand = TRUE, select="all")

# hits is a Hits object; pull out what you need
q     <- queryHits(hits)      # row in grA
s     <- subjectHits(hits)    # row in grB
d     <- mcols(hits)$distance # distance in bp (0 if they overlap)

result = cbind(as(reads[q], "DataFrame"), as(annotations[s], "DataFrame"), distance = d)
capture.output(print(result), file = files$outfile)