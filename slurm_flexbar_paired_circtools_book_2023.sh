#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 
#SBATCH --mem=8G
#SBATCH -J "flexbar"

# check if we have 3 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Read 1 file] [Read 2 file] [target dir e.g. /tmp/] [R1 marker, e.g. _R1]"
  exit
fi

# $1 -> Read 1
# $2 -> Read 2
# $3 -> Target directory


# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${1/$4/} : '\(.*\)\..*\.'`

# run on 10 CPUs
# compress with bz2
# only 30nt or longer
# no uncalled bases
# quality min phred 28
# use sanger quality values (i.e. Illumina 1.9+ encoding)
flexbar -r $1 -p $2 -t $3/$target  -n 8 -z GZ -m 30 -u 0  -q TAIL -qt 28 -as "AGATCGGAAGAG" -qf sanger -j
