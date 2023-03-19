#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -J "circtools reconstruct"

# check if we have 5 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [Sample name] [target dir e.g. /path/to/data/] [BED file] [circtools directory] [CircRNACount directory]"
  exit
fi

main_out=$2/
sample_name=$1
bed_file=$3
circtools_dir=$4
tmp_folder=/tmp/$sample_name/
circtools_out_dir=$5

mkdir -p $main_out
mkdir -p $tmp_folder

#######################

main_bam=$circtools_dir/${sample_name}.bam
main_junction=$circtools_dir/${sample_name}.Chimeric.out.junction

mate1_bam=$circtools_dir/${sample_name}.mate1.bam
mate1_junction=$circtools_dir/${sample_name}.mate1.Chimeric.out.junction

mate2_bam=$circtools_dir/${sample_name}.mate2.bam
mate2_junction=$circtools_dir/${sample_name}.mate2.Chimeric.out.junction.fixed

merged_bam=$main_out/${sample_name}_merged.bam

###################

# preprocessing

# merge both mate BAM files into one new BAM file
samtools merge -l 9 -@ 16 $merged_bam $main_bam $mate1_bam $mate2_bam

# re-index the newly aggregated BAM file
samtools index $merged_bam

circtools reconstruct -N $sample_name -D $circtools_out_dir/CircRNACount -B $merged_bam -A $bed_file -O $main_out -F $mate2_junction -R $mate2_junction -J $main_junction -T $tmp_folder -p ensembl -r 2 -e 3 -q 2 -P 16

rm /scratch/global_tmp/$sample_name/ -rf
