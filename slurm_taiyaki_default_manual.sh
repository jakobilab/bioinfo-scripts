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

module unload cuda
module load taiyaki

# check if we have 4 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [Read directory (FAST5)] [target dir e.g. /tmp/] [pretrained model] [FASTA read training set]"
  exit
fi

reads=$1
out=$2
pretrained_model=$3
fasta_file=$4
ID=$5
#CUDA_VISIBLE_DEVICES=0

# create the target directory
mkdir $out -pv

#echo "======== generate per read params"

#time generate_per_read_params.py --jobs 30 $reads > $out/modbase.tsv

#echo "======== prepare mapped reads"

#time prepare_mapped_reads.py --jobs 30 --overwrite --alphabet ACGT --mod Y A 6mA $reads $out/modbase.tsv $out/modbase.hdf5 $pretrained_model $fasta_file

#echo "======== training1"

#time train_flipflop.py --min_sub_batch_size 48 --full_filter_status --overwrite --device 0 --chunk_len_min 10000 --chunk_len_max 20000 --size 256 --stride 10 --winlen 31 --mod_factor 0.01 --outdir $out/training /biosw/taiyaki/5.0.1/models/mGru_cat_mod_flipflop.py $out/modbase.hdf5

#echo "======== training2"

#time train_flipflop.py --min_sub_batch_size 48 --full_filter_status --overwrite --device 0 --chunk_len_min 10000 --chunk_len_max 20000 --size 256 --stride 10 --winlen 31 --mod_factor 0.1 --outdir $out/training2 $out/training/model_final.checkpoint $out/modbase.hdf5

echo "======== training1"

time train_flipflop.py --full_filter_status --overwrite --device $ID --chunk_len_min 2000 --chunk_len_max 4000 --size 256 --stride 10 --winlen 31 --mod_factor 0.01 --outdir $out/default_training /biosw/taiyaki/5.0.1/models/mGru_cat_mod_flipflop.py $out/modbase.hdf5

echo "======== training2"

time train_flipflop.py --full_filter_status --overwrite --device $ID --chunk_len_min 2000 --chunk_len_max 4000 --size 256 --stride 10 --winlen 31 --mod_factor 0.1 --outdir $out/default_training2 $out/default_training/model_final.checkpoint $out/modbase.hdf5


#echo "======== basecall"
#time basecall.py --alphabet ACGT --device 0 --modified_base_output $out/basecalls.hdf5 $reads $out/training2/model_final.checkpoint  > $out/basecalls.fa

