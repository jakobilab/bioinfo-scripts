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
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 5 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [Sample name]"
  exit
fi

module load parallel

cd ${1}

#time parallel cat /scratch/tjakobi/header '>' {}_${3}.sam  :::: barcodes.txt

#time parallel time samtools view ${2} \| fgrep -w -f barcodes/barcode_{}_reads_${3}.txt '>>' {}_${3}.sam :::: barcodes.txt

#ime parallel samtools view -S -b {}_${3}.sam '>' {}_${3}.bam :::: barcodes.txt

time parallel samtools merge {}_merged.bam {}_FC1.bam  {}_FC2.bam :::: barcodes.txt

#time parallel rm {}.sam :::: barcodes.txt
