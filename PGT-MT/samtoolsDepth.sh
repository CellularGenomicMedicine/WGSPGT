###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: read depth extraction per position of the mitochondrial genome

# input: aligned WGS data from an embryo TE biopsy, chrM .bed file

# output:  read depth per position of the mitochondrial genome (.txt)

###########################################################################################################################

#!/bin/bash
#SBATCH -J MTdepth_sampleID
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5000M
#SBATCH --output=log.out
#SBATCH --error=log.err
#SBATCH --partition=research

cd dataDir

module load samtools/1.15.1

samtools depth -a -b chrM.bed -H sampleID.bam > sampleID_MTcovg.txt
