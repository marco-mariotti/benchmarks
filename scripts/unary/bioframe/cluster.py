import bioframe as bf
import pandas as pd

from lib.helpers import get_file, write_result

df = get_file("bioframe")
print(df)

df = bf.cluster(df, min_dist=None)

write_result("unary", str(df))