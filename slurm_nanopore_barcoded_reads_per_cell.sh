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
if [ ! $# == 1 ]; then
  echo "Usage: $0 [Sample name]"
  exit
fi

module load parallel

cd ${1}

mkdir barcodes

zcat real.label.gz |  cut -f2 | sort -u > barcodes.txt

#time parallel zgrep {} FC1.label.gz \| cut -f1 '>' barcodes/barcode_{}_reads_FC1.txt :::: barcodes.txt

time parallel zgrep {} real.label.gz \| cut -f1 '>' barcodes/barcode_{}_reads.txt :::: barcodes.txt


