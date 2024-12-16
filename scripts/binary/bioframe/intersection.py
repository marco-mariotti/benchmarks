import bioframe as bf
import pandas as pd


def operation(
    *,
    annotation,
    reads,
):
    df: pd.DataFrame = bf.overlap(annotation, reads)
    df = df.dropna(subset=["start_", "end_", "chrom_"])
    difference_mask = df["strand"] != df["strand_"]
    return df.loc[difference_mask]
    # if strand_behavior in {"auto", "same"}:
    #     return df.loc[difference_mask]
    # elif strand_behavior == "opposite":
    #     return df.loc[~difference_mask]
    # return df
