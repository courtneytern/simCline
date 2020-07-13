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
