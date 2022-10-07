#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Wednesday, May 4, 2016 11:14 AM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Friday, May 6, 2016 4:17 PM
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -J "pipline test"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 5 arguments
if [ ! $# == 6 ]; then
  echo "Usage: $0 [Genome Bowtie2 index argument] [Read 1 file] [Read 2 file] [target dir e.g. /awesome/project/] [Target name] {Genome Fasta}"
  exit
fi

# $1 -> rRNA index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory
# $5 -> targetname

target=$5

module load bowtie2
module load samtools

mkdir -pv $4/$target

#bowtie2 -p 20 --very-sensitive --score-min=C,-15,0 --mm -x $1 -q -1 $2 -2 $3 2> $4/$target/$target.log | samtools view -hbuS - | samtools sort -o $4/$target/${target}.bam 

#samtools view -hf 4 $4/$target/${target}.bam | samtools view -Sb - > $4/$target/${target}_unmapped.bam 

#unmapped2anchors.py $4/$target/${target}_unmapped.bam | gzip >  $4/$target/${target}_anchors.fastq.gz

bowtie2 -p 20 --score-min=C,-15,0 --reorder --mm -q -U  $4/$target/${target}_anchors.fastq.gz -x $1 | find_circ.py --genome=$6 --output=$4/$target/ --profile --throughput --name=$target > $4/$target/${target}_spliced.bed 
