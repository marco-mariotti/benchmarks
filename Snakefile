from pathlib import Path

import pandas as pd

import json
import sys

from lib.helpers import get_run_cmd, library_to_language


TIME_COMMAND = """gtime -f '{{"Elapsed": "%e", "User": "%U", "System": "%S", "MaxRSS_kB": "%M"}}' -o {output.benchmark}"""


sys.path.insert(0, os.path.abspath("lib"))

CONFIGURATION_PATH = Path("config.json")

configfile: CONFIGURATION_PATH

sample_sheet = pd.read_csv(config["sample_sheet"])

WORKDIR = Path(config["output_dir"]).expanduser()
RESULTS_DIR = Path(WORKDIR / "results")
BENCHMARK_DIR = Path(WORKDIR / "results")
DOWNLOAD_DIR = Path(WORKDIR / "downloads")

rule_files = Path("rules").rglob("**/*.smk")

for rule_file in rule_files:
    include: rule_file

annotation_files = [*sample_sheet[sample_sheet.Type == "Annotation"].Name]
read_files = [*sample_sheet[sample_sheet.Type == "Reads"].Name]

library_to_operations_mapping = pd.read_csv(config["library_to_operations_mapping"])

files_to_create = []
binary_operations = library_to_operations_mapping[library_to_operations_mapping.Type == "binary"]
for library, operation in zip(binary_operations.Library, binary_operations.Operation):
    files_to_create.extend(
        expand(
                "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/result.txt",
                RESULTS_DIR=RESULTS_DIR,
                annotation=annotation_files,
                reads=read_files,
                operation=operation,
                library=library,
            )
    )

unary_operations = library_to_operations_mapping[library_to_operations_mapping.Type == "unary"]
for library, operation in zip(unary_operations.Library, unary_operations.Operation):
    files_to_create.extend(
        expand(
                "{RESULTS_DIR}/results/unary/{operation}/{infile}/{library}/result.txt",
                RESULTS_DIR=RESULTS_DIR,
                infile=annotation_files + read_files,
                operation=operation,
                library=library,
            )
    )


rule all:
    input:
        files_to_create



def remove_suffix(p: Path) -> Path:
    return p.with_suffix("")


def get_file(name: str) -> Path:
    return DOWNLOAD_DIR / Path(sample_sheet[sample_sheet.Name == name].OutPath.iloc[0]).with_suffix("")