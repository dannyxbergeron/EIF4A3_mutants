rule rsem_reference:
    input:
        gtf = config["test_path"]['annotation'],
        fasta = config["test_path"]["genome"],
    output:
        rsem_idx = "data/references/rsem.transcripts.fa",
    params:
        base = config['test_path']['rsem_index']
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
        ref_dir = config['test_path']['rsem_index'],
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


rule agg_rsem_res:
    input:
        gene_quant = expand("results/RSEM/{id}/rsem.genes.results",
                            id=simple_id),
        trans_quant = expand("results/RSEM/{id}/rsem.isoforms.results",
                            id=simple_id)
    output:
        tpm = "results/RSEM/tpm.tsv",
        est_counts = "results/RSEM/est_counts.tsv",
        transcript_tpm = "results/RSEM/transcript_tpm.tsv",
        transcript_est_counts = "results/RSEM/transcript_est_counts.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/agg_rsem_res.py"
