

def remove_suffix(p: Path) -> Path:
    return p.with_suffix("")


def get_file(name: str) -> Path:
    return DOWNLOAD_DIR / Path(sample_sheet[sample_sheet.Name == name].OutPath.iloc[0]).with_suffix("")


rule binary:
    input:
        annotation = lambda x: get_file(x.annotation),
        reads = lambda x: get_file(x.reads),
    output:
        "{RESULTS_DIR}/results/binary/intersection/{annotation}/{reads}/{library}/{operation}.bed"
    run:
        import importlib
        module = importlib.import_module(f"scripts.binary.{wildcards.library}.{wildcards.operation}")

        benchmarked_operation = benchmark(module.operation)

        import pyranges as pr

        ann = pr.example_data.ncbi_gff

        result = benchmarked_operation(
            annotation=ann,
            reads=ann,
        )

        print(result["ResultObject"])
