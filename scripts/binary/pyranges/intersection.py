from lib.helpers import get_files, write_result

annotations, reads = get_files("pyranges")

df = annotations.overlap(reads, strand_behavior="ignore", multiple="all")

write_result("binary", df)





