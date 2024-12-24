
source("lib/helpers.R")
source("scripts/reading/genomicranges.R")

files <- get_binary_arguments()

annotations <- read_genomic_file(files$annotation_file)
reads <- read_genomic_file(files$read_file)

hits <- findOverlaps(annotations, reads, ignore.strand = TRUE)

# Get the left side with duplicates if multiple overlaps occur
# queryHits(hits) returns the indices of gr1 that overlap gr2
result <- annotations[queryHits(hits)]


capture.output(print(result), file = files$outfile)
