###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: calculate quality metrics of aligned sequencing data

# input: aligned WGS data (.bam or .cram)

# output:  sample quality metrics

###########################################################################################################################

#!/bin/bash
#SBATCH -J qualimap-sampleID
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=log.out
#SBATCH --error=log.err
#SBATCH --partition=research

module load qualimap/2.2.1

qualimap bamqc -bam sampleID.bam -nt 5
