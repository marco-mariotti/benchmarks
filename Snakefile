from pathlib import Path

import pandas as pd

import json

CONFIGURATION_PATH = Path("config.json")

CONFIG = json.load(CONFIGURATION_PATH.open())

sample_sheet = pd.read_csv(CONFIG["sample_sheet"])

WORKDIR = Path(CONFIG["output_dir"]).expanduser()
RESULTS_DIR = Path(WORKDIR / "results")
DOWNLOAD_DIR = Path(WORKDIR / "downloads")

rule_files = Path("rules").rglob("**/*.smk")

for rule_file in rule_files:
    include: rule_file





annotation_files = [*sample_sheet[sample_sheet.Type == "Annotation"].Name]
read_files = [*sample_sheet[sample_sheet.Type == "Reads"].Name]

print(
        expand(
            "{RESULTS_DIR}/results/binary/intersection/{annotation}/{reads}/intersection.bed",
            RESULTS_DIR=RESULTS_DIR,
            annotation=annotation_files,
            reads=read_files,
        )
)

rule all:
    input:
        expand(
            "{RESULTS_DIR}/results/binary/intersection/{annotation}/{reads}/intersection.bed",
            RESULTS_DIR=RESULTS_DIR,
            annotation=annotation_files,
            reads=read_files,
        )
        # Path(WORKDIR, "gencode/annotation.gtf"),
        # Path(WORKDIR, "remc/H1_cell_line.bed")

