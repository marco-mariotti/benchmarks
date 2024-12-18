import importlib
import json
from pathlib import Path
from tempfile import NamedTemporaryFile

from snakemake.script import Snakemake
from snakemake.shell import shell


from lib.benchmark import benchmark
from lib.helpers import library_to_language

snakemake: Snakemake


def run_binary_operation(
    library: str,
    operation: str,
    annotation_file: Path,
    reads_file: Path,
):
    match language := library_to_language(library=library):
        case "py":
            result, benchmarks = _run_binary_operation_python(
                library=library,
                operation=operation,
                annotation_file=annotation_file,
                reads_file=reads_file,
            )
        case "sh":
            result, benchmarks = _run_binary_operation_shell(
                library=library,
                operation=operation,
                annotation_file=annotation_file,
                reads_file=reads_file,
            )
        case "R":
            result, benchmarks = _run_binary_operation_r(
                library=library,
                operation=operation,
                annotation_file=annotation_file,
                reads_file=reads_file,
            )
        case _:
            raise NotImplementedError(
                f"There is no runner for {language} binary operations yet."
            )
    benchmarks["operation"] = operation
    benchmarks["library"] = library
    return result, benchmarks


def _run_binary_operation_python(
    library: str,
    operation: str,
    annotation_file: Path,
    reads_file: Path,
):
    operation_module = importlib.import_module(f"scripts.binary.{library}.{operation}")
    reader_module = importlib.import_module(f"scripts.reading.{library}")

    benchmarked_operation = benchmark(operation_module.operation)

    ann = reader_module.read(annotation_file)
    reads = reader_module.read(reads_file)

    result, benchmarks = benchmarked_operation(
        annotation=ann,
        reads=reads,
    )

    return result, benchmarks


def _run_binary_operation_shell(
    library: str,
    operation: str,
    annotation_file: Path,
    reads_file: Path,
):
    script_template = Path("scripts", "binary", library, f"{operation}.sh").read_text()

    benchmarked_shell = benchmark(shell)
    with NamedTemporaryFile("w+") as f:
        output_command = f" | wc -l > {f.name}"
        script = script_template.format(
            input_file1=annotation_file,
            input_file2=reads_file,
            output_command=output_command,
        )
        benchmark_results = benchmarked_shell(script)[1]
        result = Path(f.name).read_text().strip()

    return result, benchmark_results


def _run_binary_operation_r(
    library: str,
    operation: str,
    annotation_file: Path,
    reads_file: Path,
):
    script_template = Path("scripts", "binary", library, f"{operation}.sh").read_text()

    benchmarked_shell = benchmark(shell)
    with NamedTemporaryFile("w+") as f:
        output_command = f" | wc -l > {f.name}"
        script = script_template.format(
            input_file1=annotation_file,
            input_file2=reads_file,
            output_command=output_command,
        )
        benchmark_results = benchmarked_shell(script)[1]
        result = Path(f.name).read_text().strip()

    return result, benchmark_results


if __name__ == "__main__":
    library_wildcard = snakemake.wildcards.library
    operation_wildcard = snakemake.wildcards.operation

    annotation_file = Path(snakemake.input.annotation)
    reads_file = Path(snakemake.input.reads)

    output_result_file = Path(snakemake.output.result)
    output_benchmark_file = Path(snakemake.output.benchmark)

    results, benchmarks = run_binary_operation(
        library=library_wildcard,
        operation=operation_wildcard,
        annotation_file=annotation_file,
        reads_file=reads_file,
    )

    output_result_file.write_text(str(results))
    output_benchmark_file.write_text(json.dumps(benchmarks, indent=4))
