
source("lib/helpers.R")
source("scripts/reading/genomicranges.R")

files <- get_binary_arguments()

annotations <- read_genomic_file(files$annotation_file)
reads <- read_genomic_file(files$read_file)


common_seqlevels <- intersect(seqlevelsInUse(reads), seqlevelsInUse(annotations))
reads <- keepSeqlevels(reads, common_seqlevels, pruning.mode = "tidy")
annotations <- keepSeqlevels(annotations, common_seqlevels, pruning.mode = "tidy")
print(reads)
print(annotations)

fragments = subtract(annotations, reads, ignore.strand = TRUE)
result     <- unlist(fragments, use.names = FALSE) 

capture.output(print(result), file = files$outfile)