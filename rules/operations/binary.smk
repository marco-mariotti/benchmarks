

def remove_suffix(p: Path) -> Path:
    return p.with_suffix("")


def get_file(name: str) -> Path:
    return DOWNLOAD_DIR / Path(sample_sheet[sample_sheet.Name == name].OutPath.iloc[0]).with_suffix("")


rule binary:
    input:
        annotation = lambda x: get_file(x.annotation),
        reads = lambda x: get_file(x.reads),
    output:
        result = "{RESULTS_DIR}/results/binary/intersection/{annotation}/{reads}/{library}/{operation}.bed",
        benchmark = "{RESULTS_DIR}/results/binary/intersection/{annotation}/{reads}/{library}/{operation}.json",
    benchmark:
        "{RESULTS_DIR}/results/binary/intersection/{annotation}/{reads}/{library}/{operation}_benchmarks.tsv"
    run:
        import importlib
        module = importlib.import_module(f"scripts.binary.{wildcards.library}.{wildcards.operation}")

        benchmarked_operation = benchmark(module.operation)

        import pyranges as pr

        print("Reading")
        ann = pr.read_gtf(input.annotation, nrows=10000)
        print("Done reading annotatino")
        reads = pr.read_bed(input.reads)
        print("Done reading reads")

        result, benchmarks = benchmarked_operation(
            annotation=ann,
            reads=reads,
        )

        Path(output.result).write_text(str(result))
        Path(output.benchmark).write_text(json.dumps(benchmarks, indent=4))
