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
#SBATCH --mem=50G
#SBATCH -J "Nanocompore"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load nanocompore

# check if we have 5 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [Guppy folder] [Transcript file fasta] [Basecalled FASTQ file] [BAM file] [target dir e.g. /tmp/]"
  exit
fi

# create the target directory
mkdir $5 -p

# index first with nanopolish index
nanopolish index -s ${1}/sequencing_summary.txt -d ${1}/workspace ${3}

# realign raw signal to the expected reference sequence
nanopolish eventalign --threads 10 --reads ${3} --bam ${4} --genome ${2} --samples --print-read-names --scale-events --samples > ${5}/eventalign_reads.tsv

# data has to be collapsed per kmer and indexed by NanopolishComp Eventalign_collapse
NanopolishComp Eventalign_collapse -i ${5}/eventalign_reads.tsv -o ${5}/eventalign_collapsed_reads.tsv
