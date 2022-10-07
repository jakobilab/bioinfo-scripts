#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=40G
#SBATCH -J "get barcode reads"
#SBATCH -p general 

# check if we have 5 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [Sample name]"
  exit
fi

module load parallel

cd ${1}

time parallel cat /scratch/tjakobi/header '>' {}.sam  :::: barcodes.txt

time parallel time samtools view ${2} \| fgrep -w -f barcodes/barcode_{}_reads_real.txt '>>' {}.sam :::: barcodes.txt

time parallel samtools view -S -b {}.sam '>' {}.bam :::: barcodes.txt

time parallel rm {}.sam :::: barcodes.txt
