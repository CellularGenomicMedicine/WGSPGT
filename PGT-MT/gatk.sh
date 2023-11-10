###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: mitochondrial genome variant calling

# input: aligned WGS data from an embryo TE biopsy

# output:  mitochondrial variants (.g.vcf)

###########################################################################################################################

#!/bin/bash
#SBATCH -J MTGATK_sampleID
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=2000M
#SBATCH --output=log.out
#SBATCH --error=log.err
#SBATCH --partition=research

cd dataDir

module load samtools/1.15.1

samtools index sampleID.bam

module load gatk/4.3.0.0

gatk --java-options "-Xmx8g" HaplotypeCaller -R ref.fna -I sampleID.bam -O sampleID.g.vcf -ERC GVCF -L chrM
