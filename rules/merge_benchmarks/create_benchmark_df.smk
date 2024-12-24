




# rule create_benchmark_df:
#     input:
#         expand(
#             "{RESULTS_DIR}/binary/{operation}/{library}/{genome}/{nrows}/{maxlength}/benchmark.json",
#             RESULTS_DIR=RESULTS_DIR,
#             operation=operation,
#             library=library,
#         )