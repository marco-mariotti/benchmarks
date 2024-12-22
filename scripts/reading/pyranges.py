from pathlib import Path
import pandas as pd
import pyranges as pr
from pyranges import PyRanges


def read(
    f: Path,
) -> PyRanges:
    df = pd.read_table(f, names=["Chromosome", "Start", "End"], header=0)
    print(df)
    return pr.PyRanges(df)
