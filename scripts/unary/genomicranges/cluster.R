
source("lib/helpers.R")
source("scripts/reading/genomicranges.R")

files <- get_unary_arguments()

annotations <- read_genomic_file(files$annotation_file)

result <- reduce(
    annotations,
    ignore.strand = TRUE,   # FALSE replicates `bedtools cluster -s`
    with.revmap   = TRUE
)

## 2. Turn revmap into one ID per original row ---------------------
cluster_vec <- integer(length(annotations))               # pre-allocate
cluster_vec[ unlist(mcols(result)$revmap) ] <-
    rep(seq_along(result), lengths(mcols(result)$revmap))

## 3. Attach the cluster IDs to the original GRanges --------------
mcols(annotations)$cluster <- cluster_vec
result <- annotations
capture.output(print(result), file = files$outfile)