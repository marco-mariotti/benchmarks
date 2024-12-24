get_binary_arguments <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    annotation_file <- args[1]
    read_file <- args[2]
    outfile <- args[3]

    return(list(annotation_file = annotation_file,
                read_file = read_file,
                outfile = outfile))
}

get_unary_arguments <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    annotation_file <- args[1]
    outfile <- args[2]

    return(list(annotation_file = annotation_file,
                outfile = outfile))
}