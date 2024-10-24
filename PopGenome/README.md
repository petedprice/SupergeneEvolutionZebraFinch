Nextflow pipeline for running population genetics statistics

## Necessary input: 
- GATK filtered vcf files (see manuscript for filtering parameters)


## File descriptions: 
commands.bash
- Example commands for running pipeline 

metadata_full.csv
- metadata including information on contigs including length and naming
  
nextflow.config
  - nextflow configeration file
  - will need modification for personal settings

PH_main.nf
- main workflow file for running analysis
  
sample_info.csv
- sample info including sample names, karyotype, sex, etc (see headers)

## Folder descriptions: 
- scripts
  Rscripts called from modules
  
- modules
  modules called from PH_main.nf
  
- R_packages
  needed R packages


