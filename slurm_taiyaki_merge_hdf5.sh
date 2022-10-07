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
#SBATCH -c 2
#SBATCH --mem=40G
#SBATCH -J "taiyaki merge"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module unload cuda
module load taiyaki


# check if we have 6 arguments
if [ ! $# == 2 ]; then
  #echo "Usage: $0 [Read directory] [target dir e.g. /tmp/] [pretrained model] [FASTA file] [CUDA DEV] [Curlcake ID] [FASTQ file] [Reference Fasta]"
  echo "Usage: $0 [target dir e.g. /tmp/] [FAST5 dir]"
  exit
fi


out=$1
in=$2


# create the target directory
mkdir $out -pv

/biosw/taiyaki/5.0.1/misc/merge_mappedsignalfiles.py $1/modbase.hdf5 $2/*.hdf5
