#!/bin/bash

set -e

base_folder=/gpfs/home/rs1521/P_furiosus

baits=("I6V1Y9" "I6V1Y9__I6V1Y9" "I6UNQ3")

for bait in "${baits[@]}"
do
    echo "${bait} - start MSA search"
    job_id=$(qsub -q hx -v INPUT="${base_folder}/${bait}_P_furiosus_pulldown.fasta",OUTPUT="${base_folder}/${bait}_msas/",OUTPUT_LOG="${base_folder}/${bait}_output_msas.log",ERROR_LOG="${base_folder}/${bait}_error_msas.log" run_msa_search.sh)

    echo "${bait} - queue alphafold job"
    qsub -q hx -W depend=afterok:$job_id -v MSA_FOLDER="${base_folder}/${bait}_msas",OUTPUT_FOLDER="${base_folder}/${bait}_predictions",NUM_MODELS="1",NUM_RECYCLES="20" run_colabfold_predictions.sh
done
