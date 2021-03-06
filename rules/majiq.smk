rule install_majiq:
    output:
        env = directory('data/majiq_env')
    shell:
        'cd data; '
        '~/python/bin/python3.6 -m venv majiq_env && '
        'cd ../;'
        'scripts/install_majiq.sh {output.env}'

def reroute_bam(wildcards):
    id = wildcards.majiq_path.split('/')[1]
    return f'results/STAR/{id}/Aligned.sortedByCoord.out.bam'

rule star_rename:
    """
    manualy filter the GFF3 for chr3 with:
    cat Homo_sapiens.GRCh38.101.gff3
    | awk '$1 == "3" ||
           $1 == "##gff-version" ||
           $1 ~ /^#!genome/ ||
           $0 == "##sequence-region   3 1 198295559"'
           > Homo_sapiens.GRCh38.101.chr3.gff3
    """
    input:
        bam = reroute_bam,
    output:
        bam = "results/STAR/{majiq_path, [\w\d_]+/[\w\d_]+}.bam",
    params:
        chr = '3'
    conda:
        "../envs/tools.yaml"
    shell:
        'samtools index {input.bam} && '
        'samtools view -b {input.bam} {params.chr} > {output.bam}'

rule index_bam:
    input:
        bam = rules.star_rename.output.bam
    output:
        index_bam = rules.star_rename.output.bam + ".bai"
    conda:
        "../envs/tools.yaml"
    shell:
        "samtools index {input.bam}"

rule majiq_build:
    input:
        env = 'data/majiq_env',
        gff3 = config['path']['ENSEMBL_GFF3_chr3'],
        conf = 'data/majiq_build_config_{cell_type}.ini',
        star_rename_index = expand("results/STAR/{majiq_path}.bam.bai",
                                    majiq_path=majiq_path)
    output:
        out_dir = directory('results/majiq/{cell_type}/majiq_build/'),
    params:
        cmd = [
            'build '
            '--disable-ir '
            '--disable-denovo '
            '--disable-denovo-ir '
        ],
        env = 'data/majiq_env'
    threads:
        32
    shell:
        "scripts/run_majiq.sh '{params.env}' "
        "'{params.cmd}"
        "-j {threads} "
        "-o {output.out_dir} "
        "-c {input.conf} "
        "{input.gff3}' "
        "&& echo 'build done'"


def get_input(wildcards):
    res = [
        'results/majiq/{}/majiq_build/{}.majiq'.format(wildcards.cell_type, x)
        for x in rna_seq_id
        if wildcards.cell_type in x
        and wildcards.group in x
    ]
    return res

rule majiq_quant:
    input:
        majiq_build_dir = 'results/majiq/{cell_type}/majiq_build/'
    output:
        psi = 'results/majiq/{cell_type}/quant/{group}.psi.tsv'
    params:
        majiq = get_input,
        cmd = [
            'psi '
        ],
        env = 'data/majiq_env'
    threads:
        32
    shell:
        "scripts/run_majiq.sh '{params.env}' "
        "'{params.cmd} "
        "-j {threads} "
        "-o results/majiq/{wildcards.cell_type}/quant "
        "-n {wildcards.group} "
        "{params.majiq}'"


#voila cmd for tsv (could only be run in with majiq_env on...)
# voila tsv splicegraph.sql psi.voila --gene_names RPL3 -f output_file.tsv

#voila cmd for viz
# voila view splicegraph.sql deltapsi.voila



# TEST ONLY !!!!!!!
rule majiq_deltapsi:
    input:
        majiq_build_dir = 'results/majiq/stable/majiq_build/',
    output:
        deltapsi = directory('results/majiq/stable/deltapsi'),
    params:
        majiq = [
            '-grp1 ',
            'results/majiq/stable/majiq_build/stable_ctrl_1.majiq ',
            'results/majiq/stable/majiq_build/stable_ctrl_2.majiq ',
            '-grp2 ',
            'results/majiq/stable/majiq_build/stable_siRNA_eIF4A3_1.majiq ',
            'results/majiq/stable/majiq_build/stable_siRNA_eIF4A3_2.majiq ',
        ],
        cmd = [
            'deltapsi '
        ],
        env = 'data/majiq_env'
    threads:
        1
    shell:
        "scripts/run_majiq.sh '{params.env}' "
        "'{params.cmd} "
        "-j {threads} "
        "-o {output.deltapsi} "
        "-n ctrl siRNA "
        "{params.majiq}'"
