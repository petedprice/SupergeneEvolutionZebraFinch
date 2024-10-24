# How to run this pipeline for identification of allele-specific expression 


## 1. Run genotype_main.nf to create raw unfiltered SNPs
### Pipeline 
Trim reads
Allign etc 
Haplotypecaller
Salmon quantification 


## 2. Filter SNPs
script: ../var_filt/pre_wasp_Variant_filter.sh

## 3. WASP pipeline wasp_main.nf
### Pipeline 


## 4. Filter variants from wasp pipeline
script: ../var_file/Variant_filter.sh



##### DELETE ######
Before this, you will need to have filtered you genotyped vcfs from the first run, and filtered for depth, quality (using ../var_filt/asp_Variant_filter.sh).
You will then filter the genotyped vcfs from this run using (../var_file/Variant_filter.sh) before running phaser using the second set of filtered vcfs (../phaser)
