#!/bin/bash

set -e

base_folder=/gpfs/home/rs1521/P_furiosus

baits=("I6V1B6" "I6V3Z3" "I6V1B6__I6V1B6" "I6V3Z3__I6V3Z3" "I6V1B6__I6V3Z3")

for bait in "${baits[@]}"
do
    echo ${bait}
    qsub -q hx -v INPUT="${base_folder}/${bait}_P_furiosus_pulldown.fasta",OUTPUT="${base_folder}/${bait}_msas/",OUTPUT_LOG="${base_folder}/${bait}_output_msas.log",ERROR_LOG="${base_folder}/${bait}_error_msas.log" run_msa_search.sh
done
