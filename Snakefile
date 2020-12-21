import os

configfile: "config.json"

original_name = list(config['datasets'].values())
simple_id = list(config['datasets'].keys())

transiant = [x for x in simple_id if 'transiant' in x]
stable = [x for x in simple_id if 'stable' in x]

# To separate transiant and stable
majiq_path = [
    '{}/{}'.format(x.split('_')[0], x)
    for x in simple_id
]


rule all:
    input:
        expand("data/qc/{id}_fastqc.html",
                        id=simple_id),
        bam = expand("results/STAR/{id}/Aligned.sortedByCoord.out.bam",
                        id=simple_id),
        bw = expand("results/genomCov/bigwig/{id}.bw",
                        id=simple_id),
        kallisto_quant = "results/kallisto/tpm.tsv",
        rsem_quant = expand("results/RSEM/{id}/rsem.genes.results",
                            id=simple_id),
        star_rename_index = expand("results/STAR/{majiq_path}.bam.bai",
                                    majiq_path=majiq_path),
        majiq_build = expand('results/majiq/{cell_type}/majiq_build/',
                              cell_type=config['cell_types']),
        majiq_output = expand('results/majiq/{cell_type}/quant/{group}.psi.tsv',
                                cell_type=config['cell_types'], group=config['conditions']),
        # deltapsi = expand('results/majiq/{cell_type}/deltapsi/{group}.deltapsi.voila',
        #                     cell_type=config['cell_types'], group=config['conditions']),


rule download_genome:
    """ Downloads the genome from Ensembl FTP servers """
    output:
        genome = config['path']['genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gzip -d {output.genome}.gz "

rule rename_files:
    input:
        fastq = expand("data/reads/{original_name}.fastq",
                       original_name=original_name)
    output:
        new_name = expand("data/reads/{id}.fastq",
                          id=simple_id)
    run:
        for id, original in config['datasets'].items():
            old = "data/reads/{}.fastq".format(original)
            new_ = "data/reads/{}.fastq".format(id)

            print(old)
            print(new_)
            os.rename(old, new_)


# include RNA_seq
include: "rules/RNA_seq.smk"

# include genomeCov for everything
include: "rules/genomCov.smk"

# include kallisto for non-CLIP
include: "rules/kallisto.smk"

# include rsem for non-CLIP
include: "rules/rsem.smk"

# include majiq analysis
include: "rules/majiq.smk"
