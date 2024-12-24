
rule unary_python:
    wildcard_constraints:
        operation_type = "unary",
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "py"]),
        operation = "|".join(unary_operations.Operation.drop_duplicates())
    input:
        str(DOWNLOAD_DIR) + "/generated/annotation/{genome}/{nrows}/{maxlength}.bed"
    output:
        result = "{RESULTS_DIR}/{operation_type}/{operation}/{library}/{genome}/{nrows}/{maxlength}/result.txt",
        benchmark = "{RESULTS_DIR}/{operation_type}/{operation}/{library}/{genome}/{nrows}/{maxlength}/benchmark.json"
    run:
        module = "scripts.unary.{wildcards.library}.{wildcards.operation}"

        cmd = f"{TIME_COMMAND} python -m {module} {input} {output.result}"

        shell(cmd)


rule unary_shell:
    wildcard_constraints:
        operation_type = "unary",
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "sh"]),
        operation = "|".join(unary_operations.Operation)
    input:
        str(DOWNLOAD_DIR) + "/generated/annotation/{genome}/{nrows}/{maxlength}.bed"
    output:
        result = "{RESULTS_DIR}/unary/{operation}/{library}/{genome}/{nrows}/{maxlength}/result.txt",
        benchmark = "{RESULTS_DIR}/unary/{operation}/{library}/{genome}/{nrows}/{maxlength}/benchmark.json"
    run:
        script = "scripts/unary/{wildcards.library}/{wildcards.operation}"
        write_output_cmd = f" | tee >(wc -l > {output.result}.tmp) | tail > {output.result}"
        cmd = f"{TIME_COMMAND} bash {script}.sh {input[0]} {write_output_cmd}"
        shell(cmd)
        shell(f"cat {output.result}.tmp >> {output.result}")


rule unary_r:
    wildcard_constraints:
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "R"]),
        operation = "|".join(unary_operations.Operation)
    input:
        str(DOWNLOAD_DIR) + "/generated/annotation/{genome}/{nrows}/{maxlength}.bed",
    output:
        result = "{RESULTS_DIR}/unary/{operation}/{library}/{genome}/{nrows}/{maxlength}/result.txt",
        benchmark = "{RESULTS_DIR}/unary/{operation}/{library}/{genome}/{nrows}/{maxlength}/benchmark.json"
    run:
        script = "scripts/unary/{wildcards.library}/{wildcards.operation}"
        cmd = f"{TIME_COMMAND} Rscript {script}.R {input[0]} {output.result}"
        shell(cmd)