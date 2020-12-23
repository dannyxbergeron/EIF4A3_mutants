import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

kal_file = snakemake.input.kallisto
rsem_file = snakemake.input.RSEM
ref_file = snakemake.input.trans_ref

out_dir = snakemake.params.out_dir


def load_ref(file):
    df = pd.read_csv(file, sep='\t')
    df = df.loc[df.gene_name == "EIF4A2"]
    return list(df.transcript_id)

def load_and_filter(file, transcripts):

    df = pd.read_csv(file, sep='\t')
    df = df.loc[df.transcript.isin(transcripts)]
    df.set_index('transcript', inplace=True)
    df = df.T
    df.index = [
        x[:-2] if 'stable' in x
        else x
        for x in df.index
    ]
    df = df.reset_index()
    df = df.groupby('index').mean().reset_index()
    df = df.set_index('index').T

    # Filter to get only stable
    df = df[[x for x in df.columns if 'stable' in x]]
    df = df.sort_index()

    return df

def get_non_null(df1, df2):
    def filter(df):
        tmp = df.copy(deep=True)
        tmp = tmp.loc[tmp.sum(axis=1) != 0]
        return set(tmp.index)

    set1 = filter(df1)
    set2 = filter(df2)

    all_set = set1
    all_set.update(set2)

    df1_ = df1.loc[df1.index.isin(all_set)].copy(deep=True)
    df2_ = df2.loc[df2.index.isin(all_set)].copy(deep=True)

    return df1_, df2_

def graph_each(kal, rsem, col, trans_id):

    x = np.arange(len(trans_id))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(figsize=(12,8))
    rects1 = ax.bar(x - width/2, kal, width, label='kallisto')
    rects2 = ax.bar(x + width/2, rsem, width, label='RSEM')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('TPM')
    ax.set_title(f'kallisto vs RSEM for transctipt quantification {col}')
    ax.set_xticks(x)
    ax.set_xticklabels(trans_id, rotation='vertical')
    ax.legend()

    fig.tight_layout()
    plt.show()

def graph(kal_df, rsem_df):

    for col in kal_df.columns:
        graph_each(kal_df[col], rsem_df[col], col, kal_df.index)


def main():

    trans_list = load_ref(ref_file)

    kal_df = load_and_filter(kal_file, trans_list)
    rsem_df = load_and_filter(rsem_file, trans_list)

    kal_df, rsem_df = get_non_null(kal_df, rsem_df)
    print(kal_df)
    print(rsem_df)

    graph(kal_df, rsem_df)



if __name__ == '__main__':
    main()
