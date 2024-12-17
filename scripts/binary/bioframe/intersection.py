import bioframe as bf
import pandas as pd


def operation(
    *,
    annotation,
    reads,
):
    df: pd.DataFrame = bf.overlap(annotation, reads)
    return df.dropna(subset=["start_", "end_", "chrom_"])