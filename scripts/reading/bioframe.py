from pathlib import Path
import bioframe as bf


def read(
    f: Path,
    nrows: int | None = None,
):
    df = bf.read_table(
        f,
        schema=f.suffix.removeprefix("."),
        nrows=nrows,
    )

    if f.suffix == ".gtf":
        df["end"] += 1

    return df
