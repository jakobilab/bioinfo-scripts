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
#SBATCH --mem=400G
#SBATCH -J "Nanocompore sampcomp"
#SBATCH -p long
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load nanocompore

# check if we have 4 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Transcript file fasta] [YAML] [target dir e.g. /tmp/] [BED file]"
  exit
fi

# create the target directory
# mkdir $4 -p

# nanocompore sampcomp --nthreads 80 --allow_warnings --log_level debug --fasta ${1} --sample_yaml ${2} -o ${3}
#nanocompore sampcomp --allow_warnings --log_level debug --fasta ${1} --sample_yaml ${2} -o ${3} --bed ${4}

/home/tjakobi/.cache/pypoetry/virtualenvs/nanocompore-NmxgQRJS-py3.6/bin/nanocompore sampcomp --nthreads 60 --allow_warnings --log_level debug --fasta ${1} --sample_yaml ${2} -o ${3} --bed ${4}

