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
#SBATCH --gres=gpu:pascal:1

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
echo -e "SLURM process:\t$SLURM_PROCID"
echo -e "SLURM GPU ID:\t$SLURM_LOCALID"
echo -e "CUDA_DEVICE ID:\t$CUDA_DEVICE"
echo -e "CUDA_VISIBLE_DEVICES:\t$CUDA_VISIBLE_DEVICES"
echo "Device list:"
echo "$(nvidia-smi --query-gpu=name,gpu_uuid --format=csv -i $CUDA_VISIBLE_DEVICES | tail -n +2)"
echo "==== End of GPU information ===="
echo ""
#################################################################################################

module unload cuda
module load taiyaki

# check if we have 6 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [READ dir] [TRAIN dir] [OUT dir]"
  exit
fi

reads=$1
out=$3
train=$2


# create the target directory
#mkdir $out -pv

#train_flipflop.py --full_filter_status --overwrite --device $CUDA_VISIBLE_DEVICES --chunk_len_min 2000 --chunk_len_max 4000 --size 256 --stride 10 --winlen 31 --min_sub_batch_size 48 --mod_factor 0.01 --outdir $out/training /biosw/taiyaki/5.0.1/models/mGru_cat_mod_flipflop.py $out/modbase.hdf5

#train_flipflop.py --full_filter_status --overwrite --device $CUDA_VISIBLE_DEVICES --chunk_len_min 2000 --chunk_len_max 4000 --size 256 --stride 10 --winlen 31 --min_sub_batch_size 48 --mod_factor 0.1 --outdir $out/training2 $out/training/model_final.checkpoint $out/modbase.hdf5

#basecall.py --alphabet ACGT --device $CUDA_VISIBLE_DEVICES --modified_base_output $out/basecalls_full.hdf5 $reads $out/training2/model_final.checkpoint  | gzip -c > $out/basecalls_full.fa.gz
#echo XXX${CUDA_VISIBLE_DEVICES}XXX

basecall.py --alphabet ACGT --device cuda:0 --modified_base_output $out/basecalls.hdf5 $reads $train/training2/model_final.checkpoint > $out/basecalls.fa

#train_flipflop.py --device cpu --mod_factor 0.01 --outdir $out/training /biosw/taiyaki/5.0.1/models/mGru_cat_mod_flipflop.py $out/modbase.hdf5

# Second round: starting from model we just trained
#train_flipflop.py --device 0 --mod_factor 0.1 --outdir $out/training2 $out/training/model_final.checkpoint modbase.hdf5

# Basecalling
#basecall.py --device $cuda_dev --modified_base_output $out/basecalls.hdf5 $reads $out/training2/model_final.checkpoint  > $out/basecalls.fa

# compress base calls
#pigz -p 40 $out/basecalls.fa

