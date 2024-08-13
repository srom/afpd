#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=64:mem=384gb

# Variables INPUT and OUTPUT must be passed as input to qsub or uncommented below.
# e.g.: qsub -q hx -v INPUT=/gpfs/home/rs1521/abc, OUTPUT=/gpfs/home/rs1521/xyz run_msa_search.sh
#
#INPUT=/gpfs/home/rs1521/antoine/Tail_AllProtsTopredictCleanedUp.fa
#OUTPUT=/gpfs/home/rs1521/antoine/msas/

# Output and error log files can optionally be set as well:
if [ -z "$OUTPUT_LOG" ]; then
    OUTPUT_LOG="/gpfs/home/rs1521/output_run_msa_search.txt"
fi

if [ -z "$ERROR_LOG" ]; then
    ERROR_LOG="/gpfs/home/rs1521/error_run_msa_search.txt"
fi

cd $PBS_O_WORKDIR

module purge
module load ColabFold/1.5.2-foss-2022a-CUDA-11.7.0

# Running custom ColabFold version
cd /gpfs/home/rs1521/ColabFold-1.5.5

echo "Run ColabFold search"
/gpfs/easybuild/prod/software/Python/3.10.4-GCCcore-11.3.0/bin/python -m colabfold.mmseqs.search \
	--db1 uniref30_2302_db  \
	--db2 colabfold_envdb_202108_db \
	--mmseqs /gpfs/home/rs1521/MMseqs2/build/bin/mmseqs \
	--threads 64 \
	${INPUT} \
	/gpfs/home/rs1521/colabfold_database/ \
	${OUTPUT} \
	> ${OUTPUT_LOG} \
	2> ${ERROR_LOG}

echo "Rename MSA files"
cd /gpfs/home/rs1521/
module purge
. load_conda.sh
cd /gpfs/home/rs1521/amp-main

python -m src.db_utils.rename_a3m \
	--a3m_folder ${OUTPUT} \
	--fasta ${INPUT} \
	--output_folder ${OUTPUT} \
	>> ${OUTPUT_LOG} \
        2>> ${ERROR_LOG}

