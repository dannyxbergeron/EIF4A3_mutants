rule trimming:
    """ Trims the FASTQ files using Trimmomatic """
    input:
        fq = "data/reads/{id}.fastq",
    output:
        fq = "data/trimmed/{id}.fastq.gz",
    params:
        options = [
            "ILLUMINACLIP:TruSeq3-SE:2:30:10", "LEADING:5",
            "TRAILING:5", "MINLEN:36"
        ]
    log:
        "logs/trimmomatic/{id}.log"
    threads:
        32
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic SE "
        "-threads {threads} "
        "-phred33 "
        "{input.fq} "
        "{output.fq} "
        "{params.options} "
        "&> {log}"


rule qc:
    """ Assess the FASTQ quality using FastQC """
    input:
        fq = "data/trimmed/{id}.fastq.gz",
    output:
        fq_out = "data/qc/{id}_fastqc.html"
    params:
        out_dir = "data/qc"
    log:
        "logs/fastqc/{id}.log"
    threads:
         32
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.fq} "
        "&> {log}"


rule star_index:
    """ Generates the genome index for STAR """
    input:
        fasta = config["path"]["genome"],
        gtf = config["path"]['annotation']
    output:
        chrNameLength = "data/references/star_index/chrNameLength.txt"
    params:
        dir = config['path']['star_index']
    log:
        "logs/STAR/index.log"
    conda:
        "../envs/star.yaml"
    threads:
        32
    shell:
        "mkdir -p {params.dir} && "
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99"
        "&> {log}"


rule star_alignReads:
    """ Generates a bam file using STAR """
    input:
        chrNameLength = "data/references/star_index/chrNameLength.txt",
        fq = "data/trimmed/{id}.fastq.gz"
    output:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    params:
        index = config['path']['star_index'],
        output_dir = "results/STAR/{id}/"
    log:
        "logs/STAR/{id}.log"
    threads:
        32
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.index} "
        "--readFilesIn {input.fq} "
        "--runThreadN {threads} "
        "--readFilesCommand zcat "
        "--outReadsUnmapped Fastx "
        "--outFilterType BySJout "
        "--outStd Log "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.output_dir} "
        "--outFilterScoreMinOverLread 0.3 "
        "--outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 "
        "--winAnchorMultimapNmax 100 "
        "--alignEndsProtrude 5 ConcordantPair "
        "&> {log}"
