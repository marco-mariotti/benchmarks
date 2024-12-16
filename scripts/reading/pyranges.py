from pathlib import Path
import pyranges as pr
from pyranges import PyRanges


def read(
    f: Path,
    nrows: int | None = None,
) -> PyRanges:
    match f.suffix:
        case ".gtf":
            return pr.read_gtf(f, nrows=nrows)
        case ".bed":
            return pr.read_bed(f, nrows=nrows)
        case _:
            raise ValueError(f"File format not supported: {f.suffix!s}")
