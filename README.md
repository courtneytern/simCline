This repository contains  scripts and metadata files relevant to
studying the genetic history of _Drosophila simulans_ in North America, as admixed
from European and African populations.

## Getting fastq files pipeline
#### From SRA:
1. Get metadata via SRA Run Table for each study to be included
   1. Search the accession # on the SRA database
   2. Click on one entry
   3. Click on Study>All Runs. This will take you to the SRA Run Selector
   4. Click on Metadata (about halfway down the page). This will download SraRunTable.txt
   5. Change .txt to .csv to open in a more readable format
2. Parse out the relevant metadata from each SraRunTable using `parseScript.sh`
   - Will output one .csv file per study
3. Concatenate all of the study .csv files together using `concatenate.sh`
   - Output to `concatenate.csv`
4. Get the fastq files via `getfastq.sh`
   - This uses the metadata from the previous steps and `fasterqdump` to pull fastq files from the SRA database

#### From ENA:
1. Download ENA browser tools
2. Get relevant metadata into a table format
   - This can be done in any way you like; it may be easier to find the project number in SRA and get the metadata via SRA steps 1-3
3. Run `getfastqENA.sh`

Now you can go on to map the fastq files to your reference genome

## Mapping pipeline
Since this dataset uses some _simulans_ data pulled from _melanogaster_ samples, along with the raw _simulans_ data, the reference genome used in this pipeline contains both the _simulans_ and _melanogaster_ references. <br>
Both individual and pooled are mapped to the reference genome all together with `mappingscript.sh`. `mappingscriptEuro.sh` is the same thing but recongifured for the metadata format of the European DEST samples.
1. AdapterRemoval trims the adapter sequences
2. BWA MEM maps the trimmed reads to the reference genome
3. samtools sorts and indexes the mapped reads
   1. also merges paired reads, if applicable
4. Picard marks and removes duplicate reads
   - Outputs to finalmap.bam files in the specified output directory

## Pooled pipeline
1. Separate out pooled file names from individual with `pooledSetup.sh`
   1. Makes the input .txt files required for mpileup and VarScan
2. Run samtools mpileup and VarScan with `makeVCF.sh`
   - Outputs `pooledData.vcf`
3. Convert vcf to gds file with `vcf2gds.R`
   - Uses SeqArray library in R
4. Run SNP table statistics and create PCA plot in `pooledStatsPlots.R`
   - Calculates read depth. alternate allele frequencies, etc
   - More nuanced steps of PCA plot can be read in `lea.R`, but everything necessary in lea.R is in pooledStatsPlots.R.
   - Outputs `pca.pdf`
