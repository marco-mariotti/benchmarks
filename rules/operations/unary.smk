
rule unary:
    input:
        lambda x: get_file(x.infile),
    output:
        result = "{RESULTS_DIR}/results/unary/{operation}/{infile}/{library}/result.txt",
        benchmark = "{RESULTS_DIR}/results/unary/{operation}/{infile}/{library}/benchmark.json",
    benchmark:
        "{RESULTS_DIR}/results/unary/intersection/{infile}/{library}/{operation}_benchmarks.tsv"
    script:
        "scripts/run_unary_op.py"