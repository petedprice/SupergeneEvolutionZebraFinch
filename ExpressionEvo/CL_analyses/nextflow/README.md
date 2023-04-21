The first script to run is genotype_main.nf
This will quantify expression for all stated samples in the metadata.csv file. 

Pipeline 
Trim 
Allign etc 
Haplotypecaller
Salmon quantification 

Variants from this pipeline are not filtered for quality, depth etc.


2nd run: WASP
After genotype_main.nf has been run, run wasp_main.nf
Before this, you will need to have filtered you genotyped vcfs from the first run, and filtered for depth, quality (using ../var_filt/asp_Variant_filter.sh).
You will then filter the genotyped vcfs from this run using (../var_file/Variant_filter.sh) before running phaser using the second set of filtered vcfs (../phaser)
