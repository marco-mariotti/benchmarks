library(tools)  # for file_ext()
library(GenomicRanges)
library(rtracklayer)


read_genomic_file <- function(filepath) {
    ext <- tolower(file_ext(filepath))  # get the file extension in lowercase
    
    if (ext == "bed") {
        # Use rtracklayer::import() with format = "BED"
        gr <- read_bed(filepath)
    } else if (ext == "gtf") {
        # Use rtracklayer::import() with format = "GTF"
        gr <- read_gtf_raw(filepath)
    } else {
        stop("Unsupported file extension: ", ext)
    }
    
    return(gr)
}

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


read_bed <- function(bed_file) {
  gr <- import(bed_file, format = "BED")
  return(gr)
}
