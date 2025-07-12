from pathlib import Path

import pandas as pd

import json
import sys

from lib.helpers import get_run_cmd, library_to_language


sys.path.insert(0, os.path.abspath("lib"))

CONFIGURATION_PATH = Path("config.json")

configfile: CONFIGURATION_PATH

sample_sheet = pd.read_csv(config["sample_sheet"])

WORKDIR = Path(config["output_dir"]).expanduser()
RESULTS_DIR = Path(WORKDIR / "results")
BENCHMARK_DIR = Path(WORKDIR / "results")
DOWNLOAD_DIR = Path(WORKDIR / "downloads")

rule_files = Path("rules").rglob("**/*.smk")

annotation_files = [*sample_sheet[sample_sheet.Type == "Annotation"].Name]
read_files = [*sample_sheet[sample_sheet.Type == "Reads"].Name]
random_files = [*sample_sheet[sample_sheet.Type == "Random"].Name]

library_to_operations_mapping = pd.read_csv(config["library_to_operations_mapping"])
unary_operations = library_to_operations_mapping[library_to_operations_mapping.Type == "unary"]
binary_operations = library_to_operations_mapping[library_to_operations_mapping.Type == "binary"]

genome_to_parameters_mapping = pd.read_csv(config["parameters_file"])

files_to_create = []
for operation_type, library, operation in zip(
    library_to_operations_mapping.Type,
    library_to_operations_mapping.Library,
    library_to_operations_mapping.Operation,
):
    for genome, number_reads, read_maxlength in zip(
            genome_to_parameters_mapping.Genome,
            genome_to_parameters_mapping.NumberReads,
            genome_to_parameters_mapping.ReadMaxLength,
        ):
        files_to_create.extend(
            expand(
                    "{RESULTS_DIR}/{operation_type}/{operation}/{library}/{genome}/{nrows}/{maxlength}/benchmark.json",
                    RESULTS_DIR=RESULTS_DIR,
                    operation=operation,
                    library=library,
                    genome=genome,
                    nrows=number_reads,
                    maxlength=read_maxlength,
                    operation_type=operation_type,
                )
        )

for rule_file in rule_files:
    include: rule_file


rule all:
    input:
        RESULTS_DIR / "collected_results.csv"



def remove_suffix(p: Path) -> Path:
    return p.with_suffix("")


def get_file(name: str) -> Path:
    return DOWNLOAD_DIR / Path(sample_sheet[sample_sheet.Name == name].OutPath.iloc[0]).with_suffix("")