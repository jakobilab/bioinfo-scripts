#!/bin/bash

# Settings to be modified:

# gres setup:
# - gpu:tesla is currently fixed as we only have GPUs as extract resource
# - Moreover we only have cards of the tesla type
# - The number after the clon specifies the number of GPUs required,
# e.g. something between 1 and 4

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 2 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [Input folder (recursive)] [Output folder] [Kit] [Flow cell type] [model] [offset] [scale]"
  exit
fi


/scratch/ont-guppy-cpu/bin/guppy_basecaller	--trim_strategy rna\
			--verbose_logs \
			--compress_fastq \
			--fast5_out \
			-r \
			-i ${1} \
			-s ${2} \
			--model_file ${5} \
			--num_callers 64 
