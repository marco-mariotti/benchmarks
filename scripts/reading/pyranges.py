from pathlib import Path
import pandas as pd
import pyranges as pr
from pyranges import PyRanges


def read(
    f: Path,
) -> PyRanges:
    import polars as pl
    df = pl.read_csv(f, has_header=False, separator="\t").to_pandas()
    df.columns = ["Chromosome", "Start", "End"]
    # df = pd.read_table(f, names=["Chromosome", "Start", "End"], header=None)
    print(df)
    return pr.PyRanges(df)
