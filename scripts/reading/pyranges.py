from pathlib import Path
import pandas as pd
import pyranges as pr
from pyranges import PyRanges


def read(
    f: Path,
    nrows: int | None = None,
) -> PyRanges:
    match f.suffix:
        case ".gtf":
            df = pd.read_table(
                    f,
                    names=[
                        "Chromosome",
                        "Source",
                        "Feature",
                        "Start",
                        "End",
                        "Score",
                        "Strand",
                        "Frame",
                        "Attributes",
                    ],
                )
            df["End"] = df.End + 1
            return pr.PyRanges(
                df
            )
        case ".bed":
            return pr.read_bed(f, nrows=nrows)
        case _:
            raise ValueError(f"File format not supported: {f.suffix!s}")
