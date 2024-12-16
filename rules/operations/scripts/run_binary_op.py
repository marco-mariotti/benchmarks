import importlib
import json
from pathlib import Path

from snakemake.script import Snakemake
from snakemake.shell import shell


from lib.benchmark import benchmark

snakemake: Snakemake

library_wildcard = snakemake.wildcards.library
operation_wildcard = snakemake.wildcards.operation

annotation_file = snakemake.input.annotation
reads_file = snakemake.input.reads

extension = snakemake.config["library_to_language"].get(library_wildcard)

match extension:
    case "py":
        operation_module = importlib.import_module(
            f"scripts.binary.{library_wildcard}.{operation_wildcard}"
        )
        reader_module = importlib.import_module(f"scripts.reading.{library_wildcard}")

        benchmarked_operation = benchmark(operation_module.operation)

        ann = reader_module.read(Path(annotation_file), nrows=10000)
        reads = reader_module.read(Path(reads_file), nrows=10000)

        result, benchmarks = benchmarked_operation(
            annotation=ann,
            reads=reads,
        )
    case "sh":
        shell("sleep 1", bench_record=snakemake)
    case _:
        raise NotImplementedError(f"There is no runner for '.{extension}' extension binary operations yet.")


Path(snakemake.output.result).write_text(str(result))
Path(snakemake.output.benchmark).write_text(json.dumps(benchmarks, indent=4))