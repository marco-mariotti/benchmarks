from pathlib import Path

import pandas as pd

import json
import sys

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
for library, operation in zip(library_to_operations_mapping.Library, library_to_operations_mapping.Operation):
    files_to_create.extend(
        expand(
                "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/intersection.bed",
                RESULTS_DIR=RESULTS_DIR,
                annotation=annotation_files,
                reads=read_files,
                operation=operation,
                library=library,
            )
    )


rule all:
    input:
        files_to_create
