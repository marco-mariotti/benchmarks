from lib.helpers import get_file, write_result
import pyranges as pr

df = get_file("pyranges")
print(df)

df = df.merge_overlaps()

write_result("unary", str(df))