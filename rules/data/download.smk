
HUMAN_GENOME = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"

H1_CELL_LINE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM409nnn/GSM409307/suppl/GSM409307_UCSD.H1.H3K4me1.LL228.bed.gz"

files = {"gencode/annotation.gtf.gz": HUMAN_GENOME, "remc/H1_cell_line.bed.gz": H1_CELL_LINE}



rule download_hg38_comprehensive:
    output: 
        "{DATA_FOLDER}/gencode/annotation.gtf.gz"
    shell: 
        "curl -O {HUMAN_GENOME} > {output[0]}"

rule download_h1_cell_line:
    output: 
        "{DATA_FOLDER}/remc/H1_cell_line.bed.gz"
    shell: 
        "curl -O {H1_CELL_LINE} > {output[0]}"


rule gunzip:
    input: 
        "{DATA_FOLDER}/{path}.gz"
    output:
        "{DATA_FOLDER}/{path}"
    wildcard_constraints:
        path = "|".join([str(Path(p).with_suffix("")) for p in sample_sheet.OutPath])
    shell:
        "gunzip {input}"

    
# rule generate_random_data: