import pyranges as pr

def operation(
    annotation: pr.PyRanges,
    reads: pr.PyRanges,
    ):
    return annotation.intersect(reads)