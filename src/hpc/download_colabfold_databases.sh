#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -e error_colabfold_download_qsub.txt

. load_conda.sh
conda activate aria2

# Use recent version of mmseqs2
export PATH=/gpfs/home/rs1521/MMseqs2/build/bin/:$PATH

./setup_colabfold_databases.sh colabfold_database/ > output_download_colabfold_databases.txt 2> error_download_colabfold_databases.txt
