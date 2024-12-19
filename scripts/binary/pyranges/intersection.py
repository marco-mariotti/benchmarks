from lib.helpers import get_files, write_result

annotations, reads = get_files("pyranges")

df = annotations.intersect(reads, strand_behavior="ignore")

write_result("binary", str(df))





