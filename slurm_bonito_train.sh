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
#SBATCH -c 20
#SBATCH -p gpu
#SBATCH --mem=120G
#SBATCH -J "bonito train multiGPU"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de
#SBATCH --gres=gpu:tesla:4



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
module unload taiyaki
module unload guppy

module load cuda/9.2
module load bonito

# check if we have 4 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [Input dir] [Output dir]"
  exit
fi

INPUT=$1
OUTPUT=$2

# create the target directory

bonito train --amp --multi-gpu --batch 200 $OUTPUT --config /biosw/bonito/0.2.3/lib/python3.6/site-packages/bonito/models/configs/quartznet5x5_ACGTY.toml --directory $INPUT
#bonito train --amp --multi-gpu --batch 256 $OUTPUT --config /biosw/bonito/0.2.3/lib/python3.6/site-packages/bonito/models/configs/quartznet5x5.toml --directory $INPUT

