#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -J "circtools detect"

# check if we have 3 arguments
if [ ! $# == 8 ]; then
  echo "Usage: $0 [Sample sheet file] [GTF file] [Genome FASTA file] [Mate 1 file] [Mate 2 file] [BAM list file] [Repeats] [target dir e.g. project/]"
  exit
fi

circtools detect @$1 -D -an $2 -A $3 -Pi -mt1 @$4 -mt2 @$5 -B @$6 -fg -Nr 2 2 -G -k -R $7 -O $8 -F -t ${8}_circtools_temp/ -L 20 -T 8

