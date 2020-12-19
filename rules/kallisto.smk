rule create_transcriptome:
    """ Uses gffread to generate a transcriptome """
    input:
        genome = config['test_path']['genome'],
        gtf = config['test_path']['annotation']
    output:
        seqs = config['test_path']['transcriptome']
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output.seqs}"


rule generate_transcriptID_geneID:
    """
    Generating a two-column text file containing the gene -> transcript
    relationship
    """
    input:
        gtf = config['test_path']['annotation']
    output:
        map = config['test_path']['gene_name']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_transcriptID_geneID.py"


rule kallisto_index:
    """ Generates the transcriptome index for Kallisto """
    input:
        transcriptome = config['test_path']['transcriptome']
    output:
        idx = config['test_path']['kallisto_index']
    params:
        kmer = "31"
    log:
        "logs/kallisto/index.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index "
        "--index={output.idx} "
        "--kmer-size={params.kmer} "
        "{input.transcriptome} "
        "&> {log}"


rule kallisto_quant:
    """ Generates counts using Kallisto pseudo-alignment """
    input:
        idx = config['test_path']['kallisto_index'],
        fq = "data/trimmed/{id}.fastq.gz"
    output:
        quant = "results/kallisto/{id}/abundance.tsv",
        h5 = "results/kallisto/{id}/abundance.h5",
    params:
        bootstrap = "50",
        outdir = "results/kallisto/{id}"
    log:
        "logs/kallisto/{id}.log"
    threads:
        32
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "--single -l 200 -s 20 "
        "{input.fq} "
        "&> {log}"


rule combine_gene_quantification:
    """
    Custom Python script to collect and format Kallisto results for further
    processing.
    """
    input:
        datasets = expand( "results/kallisto/{id}/abundance.tsv",
                            id=simple_id),
        map = config['test_path']['gene_name']
    output:
        tpm = "results/kallisto/tpm.tsv",
        est_counts = "results/kallisto/est_counts.tsv",
        transcript_tpm = "results/kallisto/transcript_tpm.tsv",
        transcript_est_counts = "results/kallisto/transcript_est_counts.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_gene_quantification.py"
