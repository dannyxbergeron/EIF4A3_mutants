import numpy as np
import pandas as pd

gene_quant_files = snakemake.input.gene_quant
trans_quant_files = snakemake.input.trans_quant

out_tpm = snakemake.output.tpm
out_est_counts = snakemake.output.est_counts
out_trans_tpm = snakemake.output.transcript_tpm
out_trans_est_counts = snakemake.output.transcript_est_counts

def process(files, index_col):

    tpm_df = pd.DataFrame()
    count_df = pd.DataFrame()
    for file in files:
        name = file.split('/')[2]
        tmp = pd.read_csv(file, sep='\t')
        if len(tpm_df) == 0:
            tpm_df = tmp[[index_col, 'TPM']]
            tpm_df = tpm_df.rename(columns={'TPM': name})
            count_df = tmp[[index_col, 'expected_count']]
            count_df = count_df.rename(columns={'expected_count': name})
        else:
            tpm_df[name] = tpm_df[index_col].map(dict(zip(tmp[index_col],
                                                          tmp['TPM'])))
            count_df[name] = count_df[index_col].map(dict(zip(tmp[index_col],
                                                          tmp['expected_count'])))
    return tpm_df, count_df

def write(df, out_file):

    df.to_csv(out_file, sep='\t', index=False)


def main():

    t_tpm, t_count = process(trans_quant_files, 'transcript_id')
    t_tpm = t_tpm.rename(columns={'transcript_id': 'transcript'})
    t_count = t_count.rename(columns={'transcript_id': 'transcript'})

    g_tpm, g_count = process(gene_quant_files, 'gene_id')
    g_tpm = g_tpm.rename(columns={'gene_id': 'gene'})
    g_count = g_count.rename(columns={'gene_id': 'gene'})


    dfs = [t_tpm, t_count, g_tpm, g_count]
    out_files = [out_trans_tpm, out_trans_est_counts,
                out_tpm, out_est_counts]
    for df, out_file in zip(dfs, out_files):
        write(df, out_file)


if __name__ == '__main__':
    main()
