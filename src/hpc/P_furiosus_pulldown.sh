#!/bin/bash

set -e

base_folder=/gpfs/home/rs1521/P_furiosus

baits=("I6V1B6" "I6V3Z3" "I6V1B6__I6V1B6" "I6V3Z3__I6V3Z3" "I6V1B6__I6V3Z3")

for bait in "${baits[@]}"
do
    echo ${bait}
    qsub -q hx -v MSA_FOLDER="${base_folder}/${bait}_msas",OUTPUT_FOLDER="${base_folder}/${bait}_predictions",NUM_MODELS="1",NUM_RECYCLES="20" run_colabfold_predictions.sh
done
