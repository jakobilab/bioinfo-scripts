#!/usr/bin/env Rscript

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# A wrapper script for the subread R package to count features in
# mapped RNA-seq data

### imports
library(Rsubread)
library(optparse)

# organize input data
option_list = list(
 	make_option(c("-b", "--folder"), type="character", default=NULL,
			help="BAM file folder", metavar="character"),

	make_option(c("-g", "--gtf"), type="character", default=NULL,
			help="Genome for feature assignment", metavar="character"),

  make_option(c("-d", "--dir"), type="character",
      default=NULL, help="Where to save the results to?", metavar="character"),

  make_option(c("-f", "--file"), type="character",
      default="./subread.RData",
			help="Where to save the results to?", metavar="character"),

  make_option(c("-T", "--tmp"), type="character",
      default="/scratch/global_tmp/",
			help="Directory for temporary files", metavar="character"),

  make_option(c("-m", "--multi"), default=FALSE,
			help="Count multi-mapping reads?"),

  make_option(c("-t", "--threads"), type="integer", default=20,
			help="How many CPU threads should be used?", metavar="number"),

	make_option(c("-p", "--pe"), default=TRUE,
			help="Should we run in paired-end mode?"),

  make_option(c("-s", "--strand"), type="integer", default=2,
			help="Subread strandness flag", metavar="number")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (gtf)", call.=FALSE)
}

if (is.null(opt$folder)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (folder)", call.=FALSE)
}


print(opt$folder)

files <- Sys.glob(paste(opt$folder, "/*.bam",sep=""))

message("Running with the folloing BAM files:")
print(files)

subread_counts<-featureCounts(files,
# annotation
annot.inbuilt="",
annot.ext=opt$gtf,
isGTFAnnotationFile=TRUE,
GTF.featureType="gene",
#GTF.attrType="gene",
GTF.attrType="gene_id",
chrAliases=NULL,
# level of summarization
useMetaFeatures=TRUE,
# overlap between reads and features
allowMultiOverlap=TRUE,
minOverlap=1,
fracOverlap=0,
largestOverlap=FALSE,
readExtension5=0,
readExtension3=0,
read2pos=NULL,
# multi-mapping reads
countMultiMappingReads=opt$multi,
# fractional counting
fraction=FALSE,
# read filtering
minMQS=0,
splitOnly=FALSE,
nonSplitOnly=FALSE,
primaryOnly=FALSE,
ignoreDup=FALSE,
# strandness
strandSpecific=opt$strand,
# exon-exon junctions
juncCounts=FALSE,
genome=NULL,
# parameters specific to paired end reads
isPairedEnd=opt$pe,
requireBothEndsMapped=FALSE,
checkFragLength=FALSE,
minFragLength=50,
maxFragLength=600,
countChimericFragments=TRUE,
autosort=TRUE,
# number of CPU threads
nthreads=opt$threads,
# miscellaneous
maxMOp=10,
tmpDir=opt$tmp)
message("Final statistics:")
print(subread_counts$stat)

message(paste("Saving results to",opt$file))
save(subread_counts, file = opt$file)

write.table(subread_counts$counts, file = paste(opt$file,".tsv", sep = "") , row.names = TRUE, na = "", col.names = TRUE, sep = "\t", quote = FALSE)

message("Done.")
