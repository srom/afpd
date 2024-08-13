#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=384gb:ngpus=1:gpu_type=A100

# Variables MSA_FOLDER, OUTPUT_FOLDER and optionally NUM_MODELS and NUM_RECYCLES
# are provided with qsub -v argument (comma separated e.g. -v A="path1",B="path2")

if [ -z "$NUM_MODELS" ]; then
    NUM_MODELS=1
fi

if [ -z "$NUM_RECYCLES" ]; then
    NUM_RECYCLES=3
fi

# Load ColabFold library
module purge
module load ColabFold/1.5.2-foss-2022a-CUDA-11.7.0

# ColabFold 1.5.5 or more recent is needed for proper multimer support
cd /gpfs/home/rs1521/ColabFold-1.5.5

# Run ColabFold
/gpfs/easybuild/prod/software/Python/3.10.4-GCCcore-11.3.0/bin/python -m colabfold.batch \
	--num-models ${NUM_MODELS} \
	--num-recycle ${NUM_RECYCLES} \
	--model-type alphafold2_multimer_v3 \
	${MSA_FOLDER} \
	${OUTPUT_FOLDER}
