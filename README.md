# README.md

This repository contains scripts and metadata files. The goal
of this research is to elucidate the genetic history of _Drosophila simulans_
in North America, as admixed from European and African populations.
<br>
The information below explains the steps of the pipeline in order.

## Pull initial metadata from repositories
Before going through the main pipeline, get the metadata for the FASTQ files
from the repository.
#### From SRA:
Get metadata via SRA Run Table for each study to be included
1. Search the accession # on the SRA database
2. Click on one entry
3. Click on Study>All Runs. This will take you to the SRA Run Selector
4. Click on Metadata (about halfway down the page). This will download SraRunTable.txt
5. Change .txt to .csv to open in a more readable format

#### From ENA:
1. Download ENA browser tools
2. Get relevant metadata into a table format
   - This can be done in any way you like; it may be easier to find the project number in SRA and get the metadata in the same way as SRA steps 1-3

## Get FASTQ files
### Folder /1.get_data

1. `parseScript.sh`: Parse out the relevant metadata from each SraRunTable
   - Will output one .csv file per study
2. `findEuroSim.sh`: Identify which of the European DEST samples have a high enough _simulans_ contamination rate to include (>= 5 individuals)
3. `concatenate.sh`: Concatenate all of the study .csv files together (parse and Euro)
   - Output to `/metadata/concatenated.csv`
4. `getfastq.sh`: Get the FASTQ files for the data in SRA
   - This uses the metadata from the previous steps and `fasterqdump` to pull fastq files from the SRA database
5. `getfastqENA.sh`: Get the FASTQ files for the data in ENA
   - This uses ENA Browser Tools / enaGroupGet

## Pooled mapping pipeline
### Folder /2.mapping_pipeline
Since this dataset contains some _simulans_ data as decontaminated from _melanogaster_ samples, the reference genome used in this part of the pipeline includes both _simulans_ and _melanogaster_ genomes. The _simulans_ genome is separated out in downstream steps.

6. `pooledSetup.sh`: Creates all the guide files necessary for mapping the pooled samples
7. `pooledMapping.sh` and `pooledMappingEuro.sh`: Maps the pooled samples to the combined reference genome. These two scripts follow the same algorithm, but the latter is reconfigured to handle the format of the European DEST samples.
   1. `AdapterRemoval` trims the adapter sequences
   2. `BWA MEM` maps the trimmed reads to the reference genome
   3. `samtools` sorts and indexes the mapped reads; merges paired reads if applicable
   4. `Picard` marks and removes duplicate reads
      - Outputs BAM files in the specified output directory
8. `makePooledVCF.sh`: Create the VCF file with all pooled samples
   - Outputs `pooledData.vcf`

## Individual mapping pipeline
### Folder /2.mapping_pipeline/individ_GATK

9. `individSetup.sh`: Creates all the guide files necessary for the individual samples GATK pipeline. Code segments should be run separately; segments labelled by which step they come before.
10. Run GATK pipeline steps (thanks to Joaquin Nunez for these scripts)
    1. `1.run_fastqc_reads`: Run FASTQC and MultiQC on individual samples
    2. `2.Trim_and_Map.sh` and `2.x.Trim_and_Map_interleaved.sh`: Trims merged and unmerged reads, then maps to the reference genome. `2.x` is the same algorithm, but adapted for interleaved files.
    3. `3.run_multi_qualimap.sh`: Runs Qualimap
    4. `4.haplotype_caller.sh`: Adds read group information, index BAM files, runs GATK HaplotypeCaller.
    5. `5.MergeVCF_GenomicsDB.sh`: Merges VCF files
    6. `6.Genotype_GenomicsDB.sh`: Genotype Calling
    7. `7.Gather_VCFs.sh`: Gather VCFs from all chromosomes 
3. Convert vcf to gds file with `vcf2gds.R`
   - Uses SeqArray library in R
4. Run SNP table statistics and create PCA plot in `pooledStatsPlots.R`
   - Calculates read depth. alternate allele frequencies, etc
   - More nuanced steps of PCA plot can be read in `lea.R`, but everything necessary in lea.R is in pooledStatsPlots.R.
   - Outputs `pca.pdf`
