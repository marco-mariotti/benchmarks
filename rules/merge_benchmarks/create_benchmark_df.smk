



rule create_benchmark_df:
    input:
        files_to_create
    output:
        "{RESULTS_DIR}/collected_results.csv"
    run:
        rowdicts = []
        for f in input:
            f = Path(f).parent / "benchmark.json"
            print(f)
            rowdicts.append(json.loads(f.read_text()))

        df = pd.DataFrame.from_records(rowdicts)
        df.to_csv(output[0], index=None)
