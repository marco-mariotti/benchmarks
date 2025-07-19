



rule create_benchmark_df:
    input:
        files_to_create
    output:
        "{RESULTS_DIR}/collected_results.csv"
    run:
        rowdicts = []
        for f in input:
            f = Path(f).parent / "benchmark.jsonl"
            print(f)
            text = json.loads(f.read_text())
            print(text)
            rowdicts.append(text)

        df = pd.DataFrame.from_records(rowdicts)
        df = df.drop(["stderr", "cmd"], axis=1)
        df.to_csv(output[0], index=None)
