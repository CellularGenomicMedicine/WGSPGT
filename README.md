# WGSPGT
Data processing and visualisation for WGS-PGT (whole genome sequencing - preimplantation genetic testing).

_Note:_ the included bash scripts are optimised for computation on a high performance computing cluster
running a SLURM job scheduler. The specified resource requirements were optimised for the sequencing data
generated for this project.

## Data processing & QC
+ [coverage metrics from qualimap](QC/qualimap.sh)

## PGT-AO (PGT for aneuploidy origins)


## PGT-SR (PGT for structural rearrangements)
Embryo trophectoderm biopsy (and parental/reference) data was processed with the following steps:
1. The data were processed as per the PGT-M processing up to and including haplarithmisis.
2. [breakpoint analysis using Manta](PGT-SR/manta.sh)
3. [Relevant breakpoint extraction](PGT-SR/extractBreakpoints.R)
4. Custom [visualisation](PGT-SR/plotSR.R) of haplarithms including breakpoint information & chromosome schematics

## PGT-MT (PGT for mitochondrial disorders)
Embryo trophectoderm biopsy WGS data was processed with the following steps:  
1. The data were processed as per the PGT-M processing up to and including the alignment step.
    (alignment was done to the Hg38 reference genome including the mitochondrial "chromosome")
2. (Samples that were deep sequenced (30-40X) were subsampled with the aforementioned procedure.
3. Mitochondrial DNA coverage calculation & visualisation
    + [Depth per position calculation](PGT-MT/samtoolsDepth.sh)
    + [Combination of the per sample data](PGT-MT/combineDepths.R)
    + [Visualisation of the depth per position](PGT-MT/circosPlot.R)
    + [Extraction and visualisation](PGT-MT/pathogenicCoverageHistogram.R) of the [pathogenic sites (derived from MITOMAP November 2023)](PGT-MT/pathMITO.csv)
4. Heteroplasmy level calculation
    + [Mitochondrial variant calling with GATK](PGT-MT/gatk.sh)
    + [Heteroplasmy level calculation](PGT-MT/heteroplasmy.R)


### References:
__MITOMAP:__
