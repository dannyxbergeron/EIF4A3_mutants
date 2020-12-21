rule rsem_reference:
    input:
        gtf = config["path"]['annotation'],
        fasta = config["path"]["genome"],
    output:
        rsem_idx = "data/references/rsem.transcripts.fa",
    params:
        base = config['path']['rsem_index']
    log:
        "logs/RSEM/ref_log.tok"
    threads:
        32
    conda:
        "../envs/rsem.yaml"
    shell:
        "rsem-prepare-reference "
        "--star "
        "-p {threads} "
        "--gtf {input.gtf} "
        "{input.fasta} "
        "{params.base}"

rule run_rsem:
    input:
        rsem_idx = "data/references/rsem.transcripts.fa",
        fq = "data/trimmed/{id}.fastq.gz",
    output:
        quant = "results/RSEM/{id}/rsem.genes.results"
    params:
        base_name = "results/RSEM/{id}/rsem",
        ref_dir = config['path']['rsem_index'],
        tmp_dir = "results/RSEM/{id}/tmp"
    log:
        "logs/RSEM/{id}.log"
    threads:
        32
    conda:
        "../envs/rsem.yaml"
    shell:
        "rsem-calculate-expression "
        "-p {threads} "
        "--star "
        "--output-genome-bam "
        "--temporary-folder {params.tmp_dir} "
        "--star-gzipped-read-file "
        "{input.fq} "
        "{params.ref_dir} "
        "{params.base_name} &> {log}"
