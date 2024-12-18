

read_gtf_raw <- function(gtf_file) {
  df <- read.table(gtf_file, sep = "\t", comment.char = "#",
                   quote = "", header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("seqname", "source", "feature", "start", "end", 
                    "score", "strand", "frame", "attributes")
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE,
                                 seqnames.field = "seqname",
                                 start.field = "start",
                                 end.field = "end",
                                 strand.field = "strand")
  return(gr)
}