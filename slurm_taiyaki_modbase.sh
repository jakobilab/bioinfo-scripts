#!/bin/bash

# Copyright (C) 2020 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH -p general
#SBATCH --mem=200G
#SBATCH -J "taiyaki modbase"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module unload cuda
module load taiyaki
module load minimap2
module load parallel


# check if we have 6 arguments
if [ ! $# == 2 ]; then
  #echo "Usage: $0 [Read directory] [target dir e.g. /tmp/] [pretrained model] [FASTA file] [CUDA DEV] [Curlcake ID] [FASTQ file] [Reference Fasta]"
  echo "Usage: $0 [Read directory] [target dir e.g. /tmp/]"
  exit
fi

reads=$1
out=$2
pretrained_model=$3
fasta_file=$4
cuda_dev=$5
curlcake_id=$6
fastq=$7
reference_fasta=$8


# create the target directory
mkdir $out -pv


# initial mapping to assign reads to curlcakes
#time minimap2 -ax splice -uf -k14 -t 10 $reference_fasta $fastq | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > $out/initial_mapping.bam

# extract reads for our current curlcake
#samtools view $out/initial_mapping.bam | grep $curlcake_id | cut -f1 > $out/curlcake_reads.csv

# build custom reference fasta with modified sequences for curlcake reads + unmod sequences for other reads
#parallel echo ">{}" '>' $out/custom_reference.fasta\; cat  :::: $out/curlcake_reads.csv


# Create Per-read Scaling Parameters
generate_per_read_params.py --jobs 40 $reads > $out/modbase.tsv

# Modify input fasta with changed bases

# Create Mapped Read File
# We use A -> Y replacements
#prepare_mapped_reads.py  --jobs 40 --mod Y A 6mA $reads $out/modbase.tsv $out/modbase.hdf5 $pretrained_model $fasta_file 
#prepare_mapped_reads.py --limit 5  --alphabet ACGU --jobs 40 --mod Y A 6mA $reads $out/modbase.tsv $out/modbase.hdf5 $pretrained_model $fasta_file 

#prepare_mapped_reads.py --jobs 40 --overwrite --alphabet ACGT --mod Y A 6mA $reads $out/modbase.tsv $reads.hdf5 $pretrained_model $fasta_file

# Train modified base model

# Having prepared the mapped read file, the train_flipflop.py script trains a flip-flop model
# Progress is displayed on the screen and written to a log file in output directory.
# Checkpoints are regularly saved and training can be restarted from a checkpoint by replacing the model description file with the checkpoint file on the command line.

# First we use the flipflop model from Taiyaki
#train_flipflop.py --device 0 --mod_factor 0.01 --outdir $out/training taiyaki/models/mGru_cat_mod_flipflop.py modbase.hdf5

# Second round: starting from model we just trained
#train_flipflop.py --device 0 --mod_factor 0.1 --outdir $out/training2 $out/training/model_final.checkpoint modbase.hdf5

# Basecalling
#basecall.py --device $cuda_dev --modified_base_output $out/basecalls.hdf5 $reads $out/training2/model_final.checkpoint  > $out/basecalls.fa

# compress base calls
#pigz -p 40 $out/basecalls.fa

