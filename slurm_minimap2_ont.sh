#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH -J "minimap2"
#SBATCH -c 40
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de


# check if we have 6 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [Reference] [FASTQ] [target file e.g. /tmp/out.bam]"
  exit
fi


module load minimap2

time minimap2 -ax map-ont -t 40 ${1} ${2} | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > ${3}

