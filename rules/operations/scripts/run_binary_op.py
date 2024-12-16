import importlib
import json
from pathlib import Path
from lib.benchmark import benchmark

library_wildcard = snakemake.wildcards.library
operation_wildcard = snakemake.wildcards.operation

operation_module = importlib.import_module(
    f"scripts.binary.{library_wildcard}.{operation_wildcard}"
)
reader_module = importlib.import_module(f"scripts.reading.{library_wildcard}")

annotation_file = snakemake.input.annotation
reads_file = snakemake.input.reads

benchmarked_operation = benchmark(operation_module.operation)

ann = reader_module.read(Path(annotation_file), nrows=10000)
reads = reader_module.read(Path(reads_file), nrows=10000)

result, benchmarks = benchmarked_operation(
    annotation=ann,
    reads=reads,
)

Path(snakemake.output.result).write_text(str(result))
Path(snakemake.output.benchmark).write_text(json.dumps(benchmarks, indent=4))
