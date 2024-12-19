import json
from pathlib import Path
import sys

import pandas as pd


JSON = str | int | float | bool | None | dict[str, "JSON"] | list["JSON"]


def get_config() -> JSON:
    return json.load(Path("config.json").open())


def library_to_language(library: str) -> str:
    print(get_config()["library_to_language"].get(library))
    return get_config()["library_to_language"].get(library)

def get_run_cmd(library: str, command: str) -> str:
    match language := library_to_language(library=library):
        case "py":
            cmd = command.replace("/", ".")
            return f"python -m {cmd}"
        case "sh":
            return f"sh {command}.sh"
        case "R":
            return "Rscript"
    msg = f"Could not find runner for language {language}"
    raise AssertionError(msg)
    

def get_files(library: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if library == "bioframe":
        from scripts.reading.bioframe import read
    elif library == "pyranges":
        from scripts.reading.pyranges import read

    annotations = read(Path(sys.argv[1]))
    reads = read(Path(sys.argv[2]))
    return annotations, reads


def write_result(operation: str, result: str) -> None:
    if operation == "binary":
        f = sys.argv[3]
    else:
        f = sys.argv[2]

    Path(f).write_text(result)
    