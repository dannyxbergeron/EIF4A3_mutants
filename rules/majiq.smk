# rule install_majiq:
#     input:
#         bam = expand("results/STAR/ALL/{id}.bam",
#                         id=simple_id)
#     output:
#         env = directory('data/majiq_env')
#     conda:
#         '../envs/majiq.yaml'
#     shell:
#         'cd data && '
#         'virtualenv --python=python3.6.9 majiq_env && '
#         'cd ../ &&'
#         'scripts/install_majiq.sh {output.env}'

def reroute_bam(wildcards):
    id = wildcards.majiq_path.split('/')[1]
    return f'results/STAR/{id}/Aligned.sortedByCoord.out.bam'

rule star_rename:
    input:
        bam = reroute_bam
    output:
        bam = "results/STAR/{majiq_path, [\w\d_]+/[\w\d_]+}.bam"
    shell:
        'cp {input.bam} {output.bam}'

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
        gff3 = config['test_path']['ENSEMBL_GFF3'],
        conf = 'data/majiq_build_config_{cell_type}.ini',
        star_rename_index = expand("results/STAR/{majiq_path}.bam.bai",
                                    majiq_path=majiq_path)
    output:
        out_dir = directory('results/majiq/{cell_type}/majiq_build/'),
    params:
        cmd = [
            'build ',
            '--disable-ir ',
            '--disable-denovo ',
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
        for x in simple_id
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
# voila tsv splicegraph.sql psi.voila --gene_name RPL3 -f output_file.tsv

#voila cmd for viz
# voila view splicegraph.sql deltapsi.voila



# TEST ONLY !!!!!!!
# rule majiq_deltapsi:
#     input:
#         majiq_build_dir = 'results/majiq/stable/majiq_build/',
#     output:
#         deltapsi = directory('results/majiq/stable/deltapsi'),
#     params:
#         majiq = [
#             '-grp1 ',
#             'results/majiq/stable/majiq_build/stable_HBR2_1.majiq ',
#             'results/majiq/stable/majiq_build/stable_HBR2_2.majiq ',
#             '-grp2 ',
#             'results/majiq/stable/majiq_build/stable_UHR2_1.majiq ',
#             'results/majiq/stable/majiq_build/stable_UHR2_2.majiq ',
#         ],
#         cmd = [
#             'deltapsi '
#         ],
#         env = 'data/majiq_env'
#     threads:
#         4
#     shell:
#         "scripts/run_majiq.sh '{params.env}' "
#         "'{params.cmd} "
#         "-j {threads} "
#         "-o {output.deltapsi} "
#         "-n HBR UHR "
#         "{params.majiq}'"
