
source("lib/helpers.R")
source("scripts/reading/genomicranges.R")

files <- get_unary_arguments()

annotations <- read_genomic_file(files$annotation_file)

result <- reduce(annotations, ignore.strand = TRUE)


capture.output(print(result), file = files$outfile)
