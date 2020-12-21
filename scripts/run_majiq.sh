#!/bin/bash

export HTSLIB_LIBRARY_DIR=/home/danx/miniconda3/envs/snakemake/lib
export HTSLIB_INCLUDE_DIR=/home/danx/miniconda3/envs/snakemake/include

source $1/bin/activate

python data/majiq_env/lib/python3.6/site-packages/majiq/run_majiq.py $2
