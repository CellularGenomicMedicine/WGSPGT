###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: identify breakpoints in aligned sequencing data

# input: aligned WGS data (.bam)

# output:  structural variant information (.vcf)

###########################################################################################################################

#!/bin/bash
#SBATCH -J manta_sampleID
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=log.out
#SBATCH --error=log.err
#SBATCH --partition=research

mkdir dataDIR/SVCalling_manta1/

python manta/1.1.0/manta-1.1.0.centos5_x86_64/bin/configManta.py --bam dataDir/sampleID.bam --referenceFasta ref.fna --runDir dataDIR/SVCalling_manta1/

python dataDIR/SVCalling_manta1/runWorkflow.py --mode=local --jobs=20 --memGb=20

# this script can be adjusted to use manta version 1.6 by adjusting the file location of the configManta.py file to the relevant version
