#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -p medium
#SBATCH --mem=10G
#SBATCH -J "nanopore_featurecount"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# module load R

# check if we have 6 arguments
#if [ ! $# == 2 ]; then
#  echo "Usage: $0 [BAM folder] [GTF] [save to] [tmp]"
#  exit
#fi

Rscript /beegfs/homes/tjakobi/work/scripts/subread_feature_counts_single_file.R $@
 
