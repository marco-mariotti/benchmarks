
print("|".join([k for k, v in config["library_to_language"].items() if v == "py"]))
print("|".join(binary_operations.Operation))

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
        module = "scripts.binary.{wildcards.library}.{wildcards.operation}"

        cmd = f"{TIME_COMMAND} python -m {module} {input.annotation} {input.bed_file} {output.result}"

        shell(cmd)


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
        script = "scripts/binary/{wildcards.library}/{wildcards.operation}"
        write_output_cmd = f" | tee >(wc -l > {output.result}.tmp) | tail > {output.result}"
        cmd = f"{TIME_COMMAND} bash {script}.sh {input.annotation} {input.bed_file} {write_output_cmd}"
        shell(cmd)
        shell(f"cat {output.result}.tmp >> {output.result}")


rule binary_r:
    wildcard_constraints:
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "R"])
    input:
        annotation = lambda x: get_file(x.annotation),
        reads = lambda x: get_file(x.reads),
    output:
        result = "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/result.txt",
        benchmark = "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/benchmark.json",
    run:
        script = "scripts/binary/{wildcards.library}/{wildcards.operation}"
        cmd = f"{TIME_COMMAND} Rscript {script}.R {input.annotation} {input.reads} {output.result}"
        shell(cmd)