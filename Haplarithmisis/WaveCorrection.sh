#!/bin/bash
#SBATCH -J WaveCorrection
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=log.out
#SBATCH --error=log.err
#SBATCH --partition=defq

module load bioinf/bedtools/2.17.0
nucBed -bed ${sampledir}/Window10000.bed -fi ${refdir}/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.pgt.fna -o ${sampledir}/windowsize_GCcontent.txt
