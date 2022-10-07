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
#SBATCH -c 10
#SBATCH -p gpu
#SBATCH --mem=40G
#SBATCH -J "taiyaki"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de
#SBATCH --gres=gpu:tesla:2


################################################################################################
echo "==== Start of GPU information ===="

CUDA_DEVICE=$(echo "$CUDA_VISIBLE_DEVICES," | cut -d',' -f $((SLURM_LOCALID + 1)) );
T_REGEX='^[0-9]$';
if ! [[ "$CUDA_DEVICE" =~ $T_REGEX ]]; then
        echo "error no reserved gpu provided"
        exit 1;
fi

# Print debug information

echo -e "SLURM job:\t$SLURM_JOBID"
#echo -e "SLURM process:\t$SLURM_PROCID"
#echo -e "SLURM GPU ID:\t$SLURM_LOCALID"
echo -e "CUDA_DEVICE ID:\t$CUDA_DEVICE"
echo -e "CUDA_VISIBLE_DEVICES:\t$CUDA_VISIBLE_DEVICES"
echo "Device list:"
echo "$(nvidia-smi --query-gpu=name,gpu_uuid --format=csv -i $CUDA_VISIBLE_DEVICES | tail -n +2)"
echo "==== End of GPU information ===="
echo ""
#################################################################################################


module unload cuda
module load taiyaki
module load minimap2
module load parallel




nvidia-smi --query-gpu=name,gpu_uuid --format=csv -i $CUDA_VISIBLE_DEVICES

echo $CUDA_DEVICE


exit(0)
reads=$1
out=$2
cuda_dev=$3

#out=$1
#cuda_dev=$2
#reads=$3


# create the target directory
mkdir $out -pv


# initial mapping to assign reads to curlcakes
#time minimap2 -ax splice -uf -k14 -t 10 $reference_fasta $fastq | samtools view -bh -F 2324 -q 10 | samtools sort -O bam > $out/initial_mapping.bam

# extract reads for our current curlcake
#samtools view $out/initial_mapping.bam | grep $curlcake_id | cut -f1 > $out/curlcake_reads.csv

# build custom reference fasta with modified sequences for curlcake reads + unmod sequences for other reads
#parallel echo ">{}" '>' $out/custom_reference.fasta\; cat  :::: $out/curlcake_reads.csv


# Create Per-read Scaling Parameters
#generate_per_read_params.py --jobs 40 $reads > $out/modbase.tsv

# Modify input fasta with changed bases

# Create Mapped Read File
# We use A -> Y replacements
#prepare_mapped_reads.py  --jobs 40 --mod Y A 6mA $reads $out/modbase.tsv $out/modbase.hdf5 $pretrained_model $fasta_file 
#prepare_mapped_reads.py --limit 5  --alphabet ACGU --jobs 40 --mod Y A 6mA $reads $out/modbase.tsv $out/modbase.hdf5 $pretrained_model $fasta_file 
#echo prepare_mapped_reads.py --overwrite --alphabet ACGU --mod Y A 6mA $reads $out/modbase.tsv $reads.hdf5 $pretrained_model $fasta_file
#prepare_mapped_reads.py --overwrite --alphabet ACGU --mod Y A 6mA $reads $out/modbase.tsv $reads.hdf5 $pretrained_model $fasta_file

# Train modified base model

# Having prepared the mapped read file, the train_flipflop.py script trains a flip-flop model
# Progress is displayed on the screen and written to a log file in output directory.
# Checkpoints are regularly saved and training can be restarted from a checkpoint by replacing the model description file with the checkpoint file on the command line.

# First we use the flipflop model from Taiyaki
#train_flipflop.py --device 0 --mod_factor 0.01 --outdir $out/training taiyaki/models/mGru_cat_mod_flipflop.py modbase.hdf5

# see https://github.com/nanoporetech/taiyaki/issues/32#issuecomment-517638751
# for parameters

#python -m torch.distributed.launch --nproc_per_node=4 /beegfs/biosw/taiyaki/5.0.1/bin/train_flipflop.py --overwrite --chunk_len_min 2000 --chunk_len_max 4000 --size 256 --stride 10 --winlen 31 --min_sub_batch_size 48 --mod_factor 0.01 --outdir $out/training /biosw/taiyaki/5.0.1/models/mGru_cat_mod_flipflop.py $out/modbase_new.hdf5 

train_flipflop.py --full_filter_status --overwrite --device $CUDA_VISIBLE_DEVICES --chunk_len_min 2000 --chunk_len_max 4000 --size 256 --stride 10 --winlen 31 --min_sub_batch_size 48 --mod_factor 0.01 --outdir $out/training /biosw/taiyaki/5.0.1/models/mGru_cat_mod_flipflop.py $out/modbase.hdf5

train_flipflop.py --full_filter_status --overwrite --device $CUDA_VISIBLE_DEVICES --chunk_len_min 2000 --chunk_len_max 4000 --size 256 --stride 10 --winlen 31 --min_sub_batch_size 48 --mod_factor 0.1 --outdir $out/training2 $out/training/model_final.checkpoint $out/modbase.hdf5

basecall.py --alphabet ACGT --device $CUDA_VISIBLE_DEVICES --modified_base_output $out/basecalls.hdf5 $reads $out/training2/model_final.checkpoint  > $out/basecalls.fa



#train_flipflop.py --device cpu --mod_factor 0.01 --outdir $out/training /biosw/taiyaki/5.0.1/models/mGru_cat_mod_flipflop.py $out/modbase.hdf5

# Second round: starting from model we just trained
#train_flipflop.py --device 0 --mod_factor 0.1 --outdir $out/training2 $out/training/model_final.checkpoint modbase.hdf5

# Basecalling
#basecall.py --device $cuda_dev --modified_base_output $out/basecalls.hdf5 $reads $out/training2/model_final.checkpoint  > $out/basecalls.fa

# compress base calls
#pigz -p 40 $out/basecalls.fa

