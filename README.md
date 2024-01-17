# WGSPGT
Data processing and visualisation for WGS-PGT (whole genome sequencing - preimplantation genetic testing).

The WGS-PGT approach, for the first time, enables all forms of PGT on nuclear and mitochondrial DNA in a single-assay.
Specifically, our innovative approach outperforms traditional and state-of-the-art PGT, including:
1. PGT for genetic indications in complex genomic regions;
2. direct detection of single- and few-base pair genetic variations
3. novel form of PGT-A that uncovers segregational origin (meiotic vs. mitotic) of aneuploidies and their level of mosaicism, an approach we coin PGT-AO
4. (in)direct detection of the translocation breakpoints and inheritance of normal and derivative chromosomes
5. PGT for mitochondrial DNA (mtDNA) disorders

_Note:_ the included bash scripts are optimised for computation on a high performance computing cluster
running a SLURM job scheduler. The specified resource requirements were optimised for the sequencing data
generated for this project.

## Haplarithmisis
### Pre-processing
Sequencing output is processed via the following steps (scripts not included here)
1. Demultiplexing
2. FASTQC
3. Alignment using BWA-MEM2
4. BAMQC
5. GATK Joint-genotyping using GATK4 <br>
   Family-wise batches (i.e. multisample vcf files) are generated using gatk-haplotype-joint calling. The resulting vcf file was analysed by Haplarithmisis. <br>
   _GATK4 from Broad Institute is available [here](https://gatk.broadinstitute.org/hc/en-us)_

   
### Running Haplarithmisis
These are all relevant R scripts to process whole genome sequencing data using Haplarithmisis for WGS-PGT 
1. [MetaInfo](Haplarithmisis/MetaInfo) <br>
   input: CSV file with family information [(example attached)](Haplarithmisis/ExampleSamplesheet.csv)
   + Sample ID (_Sample ID of each of the family members / embryos_)
   + Family number (_PGT + familynumber_)
   + PGT (_diagnostics / research_)
   + Sample Status (_E = Embryo, Combination of U = unaffected or A = affected and family member: F = Father, M = Mother, S = Sibling, GF = GrandFather, GM = GrandMother_)
   + Family interval (_chr_startposition_endposition_parent_, _example: chr2_1001_1002_Pat_)  
   + Family second interval (_in case of a compound heterozygous mutation, if not applicable: Non Defined (ND)_)
   + Family indication (_GENE_ + "_ _PGT"_)
   + Family Dnr (_prefix D"year of analysis/number"_)
   
   input: PGT config file (.txt) with path to scripts, samplesheet and default parameters
   + default parameters: Win=10, gammaBAF=10, bin=10000, Window=22, gammaSC=300, gammaMC=50, ExtInt=1, plateau=100, gtypemodulator_window=10000
   
3. [ConvertGenotype](Haplarithmisis/ConvertGenotype)
4. [QDNASeq](Haplarithmisis/QDNASeq)
#### EmbryoTest: when embryo sequencing information is present (continue with step 5 NucBedPrep)
5. [NucBedPrep](Haplarithmisis/NucBedPrep)
6. [PGT Wave correction](Haplarithmisis/WaveCorrection.sh)
7. [Haplarithmisis](Haplarithmisis/Haplarithmisis)
8. [EmbryoTestReportData](Haplarithmisis/EmbryoTestReportData)
9. [EmbryoTestReportPlot](Haplarithmisis/EmbryoTestReportPlot)

#### PreTest: if no embryo sequencing information is present (continue with step 5 PreTestReportData)
5. [PreTestReportData](Haplarithmisis/PreTestReportData)
6. [PreTestReportPlot](Haplarithmisis/PreTestReportPlot)

+ [functions](Haplarithmisis/functions)
+ [Rda](Haplarithmisis/Rda)

## Data processing & QC
+ [Subsampling](QC/SubSampling.sh) to desired target coverage
+ [Visualisation of lab step timings](QC/labtimings.R)
+ [coverage metrics from qualimap](QC/qualimap.sh)
+ [Extract qualimap output](QC/QualimapExtraction.R)
+ [Visualisation coverage metrics](QC/CoveragePlots.R)
+ Visualisation of Mendelian inconsistency for [validation](QC/AutosomalMendIncValidationPlot.R) and [pilot per subsampled target coverage and validation per chromosome](QC/MendInconsistencyPlotsSupplement.R)
+ Visualisation of Haplotype concordance for [pilot at subsampled target coverages](QC/HaplotypeConcordancePilotPlot.R) and [validation](QC/HaplotypeConcordanceValidationPlot.R)
+ [liftover coordinates from onePGT output](QC/LiftOverhg19.R)
+ informative SNP [binning](QC/binning.R) and [chromosome heatmap visualisation](QC/coverageHeatmap.R) - see PGT-SR folder for chromosome coordinate scripts.

## PGT-M (PGT for monogenic disorders)
5. [Haplarithm plotting](PGT-M/haplarithm.R)
6. [direct mutation analysis visualisation](PGT-M/directMutationPlot.R)


## PGT-AO (PGT for aneuploidy origins)
1. input: CSV file with family information [(example attached for parents-only haplarithmisis)](Haplarithmisis/ExampleSamplesheet_parentsOnlyPGTA.csv)
2. [Haplarithm plotting](PGT-AO/haplarithm.R) with [chromosome ideogram](PGT-AO/ideogram.R) - see PGT-SR folder for ideogram coordinate scripts.


## PGT-SR (PGT for structural rearrangements)
Embryo trophectoderm biopsy (and parental/reference) data was processed with the following steps:
1. The data were processed as per the PGT-M processing up to and including haplarithmisis.
2. Deep (30-40X) sequenced data was subsampled as per PGT-M and the (segmented) logRs were [plotted](PGT-SR/plotLogRs.R)
3. [breakpoint analysis using Manta](PGT-SR/manta.sh)
4. [Relevant breakpoint extraction](PGT-SR/extractBreakpoints.R)
5. Custom [visualisation](PGT-SR/plotSR.R) of haplarithms including breakpoint information & chromosome schematics
       + generation of chromosome fill / outline coordinates for [normal](PGT-SR/normalCoordinates.R) and [affected](PGT-SR/affectedCoordinates.R)
6. Visualisation of [copy number variation](PGT-SR/Veriseq_ggplot.R) from VeriSeq output.

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


#### References:
__MITOMAP:__ https://www.mitomap.org/foswiki/bin/view/MITOMAP/ConfirmedMutations
