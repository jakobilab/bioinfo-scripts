#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=5G
#SBATCH -J "get barcode reads"
#SBATCH -p medium 

# check if we have 5 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [Sample name]"
  exit
fi


cat /scratch/tjakobi/header > ${1}.sam

time samtools view ${2} | fgrep -w -f barcodes/barcode_${1}_reads_real.txt >> ${1}.sam

samtools view -S -b ${1}.sam > ${1}.bam

time parallel rm {}.sam :::: barcodes.txt
