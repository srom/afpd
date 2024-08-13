#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=64gb:ngpus=1:gpu_type=A100
#PBS -e error_run_colabfold_batch_antoine.txt
#PBS -o output_run_colabfold_batch_antoine.txt

cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

cat $PBS_NODEFILE

module purge
module load ColabFold/1.5.2-foss-2022a-CUDA-11.7.0

colabfold_batch \
	--num-models 3 \
	--num-recycle 5 \
	--templates \
	antoine/clusters_rep_seq_70_300 \
	antoine/outputs/

