

rule binary:
    input:
        annotation = lambda x: get_file(x.annotation),
        reads = lambda x: get_file(x.reads),
    output:
        result = "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/result.txt",
        benchmark = "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/benchmark.json",
    benchmark:
        "{RESULTS_DIR}/results/binary/intersection/{annotation}/{reads}/{library}/{operation}_benchmarks.tsv"
    script:
        "scripts/run_binary_op.py"