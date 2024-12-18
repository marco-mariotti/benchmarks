import bioframe as bf
import pandas as pd


def operation(
    *,
    df: pd.DataFrame,
):
    return bf.merge(df)
