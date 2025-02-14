import bioframe as bf
import pandas as pd

from lib.helpers import get_files, write_result

annotations, reads = get_files("bioframe")

df: pd.DataFrame = bf.overlap(annotations, reads)
df = df.dropna(subset=["start_", "end_", "chrom_"]).drop(
    columns=[c for c in df.columns if c.endswith("_")], axis="columns"
)

write_result("binary", str(len(df)))
