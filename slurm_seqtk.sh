#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -J "seqtk"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load seqtk

# check if we have 2 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [Read file] [Percent] [target dir e.g. /tmp/]"
  exit
fi

# $1 -> Read
# $2 -> Percent, i.e. 0.01
# $3 -> Target directory

time seqtk sample -s10 $1 $2 | gzip > $3/$1 

