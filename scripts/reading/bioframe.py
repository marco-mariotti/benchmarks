from pathlib import Path
import pandas as pd
import bioframe as bf


def read(
    f: Path,
    nrows: int | None = None,
):
    df = bf.read_table(
        f,
        schema="bed3",
        nrows=nrows,
    )

    if f.suffix == ".gtf":
        df["end"] += 1

    return df
