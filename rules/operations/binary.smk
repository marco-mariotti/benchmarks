

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
    script:
        "scripts/run_binary_op.py"