Nextflow pipeline for running population genetics statistics

## Pre workflow to do 
Filter VCFs using GATK (see manuscript for description) 
Generate degeneracy-all-sites.bed using degenotate from https://github.com/harvardinformatics/degenotate to run stats on four-fold (FD) and zero-fold (0D) degenerate sites. 


## Necessary input: 
nextflow run $projectdir/SupergeneEvolutionZebraFinch/PopGenome/PH_main.nf \
	--metadata $projectdir/SupergeneEvolutionZebraFinch/PopGenome/metadata_full.csv \   
	--degen $datadir/degeneracy-all-sites.bed \
	--vcfs $datadir/vcfs/ \
	--sample_info $projectdir/SupergeneEvolutionZebraFinch/PopGenome/sample_info.csv \
	-resume



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

## Workflow description 
ZDed=ZD_BED(contig_ch)
	FDed=FD_BED(contig_ch)
        Xged=Xgene_BED(contig_ch)
	NCed=NCgene_BED(contig_ch)
        NCXGed=NCXG_BED(NCed.combine(Xged, by: [0,1,2]))

	VCF_filtered=vcftools_filt(contig_ch)

        pseudohaploided=pseudohaploid(VCF_filtered)

	VCF_nobeded=tabix(pseudohaploided)

	VCF_bededZ=vcftools_bedZD(ZDed.combine(pseudohaploided, by: [0,1,2]))
        VCF_bededD=vcftools_bedFD(FDed.combine(pseudohaploided, by: [0,1,2]))
        VCF_bededxG=vcftools_bedxG(Xged.combine(pseudohaploided, by: [0,1,2]))
        VCF_bededNC=vcftools_bedNC(NCed.combine(pseudohaploided, by: [0,1,2]))
        VCF_bededNCXG=vcftools_bedNCXG(NCXGed.combine(pseudohaploided, by: [0,1,2]))


	PopGenomedZ=PopGenomeZD(VCF_bededZ)
        PopGenomedF=PopGenomeFD(VCF_bededD)
        PopGenomedxG=PopGenomexG(VCF_bededxG)
        PopGenomedNCXG=PopGenomeNCXG(VCF_bededNCXG)
	PopGenomeAG=PopGenomeAG(VCF_nobeded)
        PopGenomeNC=PopGenomeNC(VCF_bededNC)

