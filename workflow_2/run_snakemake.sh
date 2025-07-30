#!/bin/bash
#SBATCH -p main
#SBATCH -n10
#SBATCH --time=24:00:00

source ~/miniconda3/etc/profile.d/conda.sh

conda activate target_search

snakemake -p all --use-conda --cores 10 --latency-wait 60