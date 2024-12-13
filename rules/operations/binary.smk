

def remove_suffix(p: Path) -> Path:
    return p.with_suffix("")


def get_file(name: str) -> Path:
    return DOWNLOAD_DIR / Path(sample_sheet[sample_sheet.Name == name].OutPath.iloc[0]).with_suffix("")


rule intersect:
    input:
        annotation = lambda x: get_file(x.annotation),
        reads = lambda x: get_file(x.reads),
    output:
        "{RESULTS_DIR}/results/binary/intersection/{annotation}/{reads}/intersection.bed"
    run:
        import pyranges
        annotation = pr.read_gtf(input.annotation)
        print(annotation)
        reads = pr.read_bed(input.reads)
