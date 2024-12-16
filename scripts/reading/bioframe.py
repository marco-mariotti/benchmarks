from pathlib import Path
import bioframe as bf


def read(
    f: Path,
    nrows: int | None = None,
):
    return bf.read_table(
        f,
        schema=f.suffix.removeprefix("."),
        nrows=nrows,
    )
