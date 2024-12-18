import pyranges as pr


def operation(
    *,
    df: pr.PyRanges,
):
    return df.merge_overlaps()
