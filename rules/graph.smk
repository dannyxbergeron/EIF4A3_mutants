rule graph_isoforms:
    input:
        kallisto = 'CLUSTER_RESULTS/kallisto/transcript_est_counts.tsv',
        RSEM = 'CLUSTER_RESULTS/RSEM/transcript_est_counts.tsv',
        trans_ref = 'data/references/transcriptsId_geneId_geneName.tsv'
    output:
        'tok'
    params:
        out_dir = 'CLUSTER_RESULTS/'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/graph_isoforms.py"
