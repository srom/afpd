#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=8:mem=384gb:ngpus=1:gpu_type=A100
#PBS -e error_colabfold_heptamer.txt
#PBS -o output_colabfold_heptamer.txt

cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

cat $PBS_NODEFILE

module purge
module load ColabFold/1.5.2-foss-2022a-CUDA-11.7.0

# Running custom ColabFold version
cd /gpfs/home/rs1521/ColabFold-1.5.5

/gpfs/easybuild/prod/software/Python/3.10.4-GCCcore-11.3.0/bin/python -m colabfold.batch \
	--num-models 3 \
	--num-recycle 20 \
	--model-type alphafold2_multimer_v3 \
	/gpfs/home/rs1521/chahrazad/msa_pgaptmp_000187_2_heptamer/ \
	/gpfs/home/rs1521/chahrazad/predictions_pgaptmp_000187_2_heptamer/

