#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 64
#SBATCH -w galatea
#SBATCH --mem=400G
#SBATCH -p long
#SBATCH -J "FAST5-convert"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 6 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [input folder] [output folder] [reads per multi file]"
  exit
fi

single_to_multi_fast5 -i ${1} -s ${2} -t 64 -n ${3} --recursive

