import importlib
import json
from pathlib import Path
from tempfile import NamedTemporaryFile

from snakemake.script import Snakemake
from snakemake.shell import shell


from lib.benchmark import benchmark
from lib.helpers import library_to_language

snakemake: Snakemake


def run_unary_operation(
    library: str,
    operation: str,
    infile: Path,
):
    match language := library_to_language(library=library):
        case "py":
            result, benchmarks = _run_unary_operation_python(
                library=library,
                operation=operation,
                infile=infile,
            )
        case "sh":
            result, benchmarks = _run_unary_operation_shell(
                library=library,
                operation=operation,
                infile=infile,
            )
        case _:
            raise NotImplementedError(
                f"There is no runner for {language} unary operations yet."
            )
    benchmarks["operation"] = operation
    benchmarks["library"] = library
    return result, benchmarks


def _run_unary_operation_python(
    library: str,
    operation: str,
    infile: Path,
):
    operation_module = importlib.import_module(f"scripts.unary.{library}.{operation}")
    reader_module = importlib.import_module(f"scripts.reading.{library}")

    benchmarked_operation = benchmark(operation_module.operation)

    df = reader_module.read(infile)

    result, benchmarks = benchmarked_operation(
        df=df,
    )

    return result, benchmarks


def _run_unary_operation_shell(library: str, operation: str, infile: Path):
    script_template = Path("scripts", "unary", library, f"{operation}.sh").read_text()

    benchmarked_shell = benchmark(shell)
    with NamedTemporaryFile("w+") as f:
        output_command = f" | wc -l > {f.name}"
        script = script_template.format(
            infile=infile,
            output_command=output_command,
        )
        benchmark_results = benchmarked_shell(script)[1]
        result = Path(f.name).read_text().strip()

    return result, benchmark_results


if __name__ == "__main__":
    library_wildcard = snakemake.wildcards.library
    operation_wildcard = snakemake.wildcards.operation

    infile = Path(snakemake.input[0])

    output_result_file = Path(snakemake.output.result)
    output_benchmark_file = Path(snakemake.output.benchmark)

    results, benchmarks = run_unary_operation(
        library=library_wildcard,
        operation=operation_wildcard,
        infile=infile,
    )

    output_result_file.write_text(str(results))
    output_benchmark_file.write_text(json.dumps(benchmarks, indent=4))
