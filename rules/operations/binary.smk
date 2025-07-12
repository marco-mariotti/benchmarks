from helpers import run_case 


print("|".join([k for k, v in config["library_to_language"].items() if v == "py"]))
print("|".join(binary_operations.Operation))

#     number_rows: int,
#     library: str,
#     genome: str,
#     max_length: int,

rule binary_python:
    input:
        annotation = str(DOWNLOAD_DIR) + "/generated/annotation/{genome}/{nrows}/{maxlength}.bed",
        bed_file = str(DOWNLOAD_DIR) + "/generated/reads/{genome}/{nrows}/{maxlength}.bed"
    wildcard_constraints:
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "py"]),
        operation = "|".join(binary_operations.Operation.drop_duplicates())
    output:
        result = "{RESULTS_DIR}/binary/{operation}/{library}/{genome}/{nrows}/{maxlength}/result.txt",
        benchmark = "{RESULTS_DIR}/binary/{operation}/{library}/{genome}/{nrows}/{maxlength}/benchmark.json"
    run:
        module = f"scripts.binary.{wildcards.library}.{wildcards.operation}"
        print(subprocess.call("which python", shell=True))
        print(subprocess.call("python --version", shell=True))

        cmd = f"python -m {module} {input.annotation} {input.bed_file} {output.result}"

        run_case(
            cmd,
            benchmark_file=output.benchmark,
            result_file=output.result,
            number_rows=wildcards.nrows,
            library=wildcards.library,
            genome=wildcards.genome,
            max_length=wildcards.maxlength,
            operation=wildcards.operation,
        )




rule binary_shell:
    wildcard_constraints:
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "sh"]),
        operation = "|".join(binary_operations.Operation)
    input:
        annotation = str(DOWNLOAD_DIR) + "/generated/annotation/{genome}/{nrows}/{maxlength}.bed",
        bed_file = str(DOWNLOAD_DIR) + "/generated/reads/{genome}/{nrows}/{maxlength}.bed"
    output:
        result = "{RESULTS_DIR}/binary/{operation}/{library}/{genome}/{nrows}/{maxlength}/result.txt",
        benchmark = "{RESULTS_DIR}/binary/{operation}/{library}/{genome}/{nrows}/{maxlength}/benchmark.json"
    run:
        script = f"scripts/binary/{wildcards.library}/{wildcards.operation}"
        write_output_cmd = f" | tail > {output.result}"
        cmd = f"bash {script}.sh {input.annotation} {input.bed_file} {write_output_cmd}"

        run_case(
            cmd,
            benchmark_file=output.benchmark,
            result_file=output.result,
            number_rows=wildcards.nrows,
            library=wildcards.library,
            genome=wildcards.genome,
            max_length=wildcards.maxlength,
            operation=wildcards.operation,
        )


rule binary_r:
    wildcard_constraints:
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "R"]),
        operation = "|".join(binary_operations.Operation)
    input:
        annotation = str(DOWNLOAD_DIR) + "/generated/annotation/{genome}/{nrows}/{maxlength}.bed",
        bed_file = str(DOWNLOAD_DIR) + "/generated/reads/{genome}/{nrows}/{maxlength}.bed"
    output:
        result = "{RESULTS_DIR}/binary/{operation}/{library}/{genome}/{nrows}/{maxlength}/result.txt",
        benchmark = "{RESULTS_DIR}/binary/{operation}/{library}/{genome}/{nrows}/{maxlength}/benchmark.json"
    run:
        script = f"scripts/binary/{wildcards.library}/{wildcards.operation}"
        cmd = f"Rscript {script}.R {input.annotation} {input.bed_file} {output.result}"
        run_case(
            cmd,
            benchmark_file=output.benchmark,
            result_file=output.result,
            number_rows=wildcards.nrows,
            library=wildcards.library,
            genome=wildcards.genome,
            max_length=wildcards.maxlength,
            operation=wildcards.operation,
        )