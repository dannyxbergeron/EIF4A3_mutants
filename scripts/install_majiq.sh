#!/bin/bash

export HTSLIB_LIBRARY_DIR=/home/berd2710/miniconda3/envs/snakemake/lib
export HTSLIB_INCLUDE_DIR=/home/berd2710/miniconda3/envs/snakemake/include

source $1/bin/activate

pip install git+https://bitbucket.org/biociphers/majiq_academic.git#egg=majiq
