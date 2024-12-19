

rule binary_python:
    wildcard_constraints:
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "py"])
    input:
        annotation = lambda x: get_file(x.annotation),
        reads = lambda x: get_file(x.reads),
    output:
        result = "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/result.txt",
        benchmark = "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/benchmark.json",
    run:
        module = "scripts.binary.{wildcards.library}.{wildcards.operation}"

        cmd = f"{TIME_COMMAND} python -m {module} {input.annotation} {input.reads} {output.result}"

        shell(cmd)


rule binary_shell:
    wildcard_constraints:
        library = "|".join([k for k, v in config["library_to_language"].items() if v == "sh"])
    input:
        annotation = lambda x: get_file(x.annotation),
        reads = lambda x: get_file(x.reads),
    output:
        result = "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/result.txt",
        benchmark = "{RESULTS_DIR}/results/binary/{operation}/{annotation}/{reads}/{library}/benchmark.json",
    run:
        script = "scripts/binary/{wildcards.library}/{wildcards.operation}"
        write_output_cmd = f" | tee >(wc -l >> {output.result}) | tail > {output.result}"
        cmd = f"{TIME_COMMAND} bash {script}.sh {input.annotation} {input.reads} {write_output_cmd}"
        shell(cmd)
    