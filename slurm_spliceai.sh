#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1 
#SBATCH --mem=50G
#SBATCH -J "spliceAI"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de
#SBATCH -p gpu

#module load spliceai

# check if we have 2 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [Input file] [Output file] [Genome (fasta)] [Annotation] [GPU]"
  exit
fi

export CUDA_VISIBLE_DEVICES=${5}

# spliceai -I /scratch/tjakobi/test.vcf -O /scratch/tjakobi/out.vcf -A grch38 -R /scratch/tjakobi/hg38.fa 

spliceai 	-I ${1} \
		-D 50 \
		-O ${2} \
		-A ${4} \
		-R ${3}

