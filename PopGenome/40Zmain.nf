nextflow.enable.dsl=2

include { FD_BED } from './modules/FD_BED.nf'
include { ZD_BED } from './modules/ZD_BED.nf'
include { Xgene_BED } from './modules/Xgene_BED.nf'
include { NCgene_BED } from './modules/NCgene_BED.nf'
include { tabix } from './modules/tabix.nf'
include { vcftools_filt } from './modules/vcftools_filt.nf'
include { vcftools_bed as vcftools_bedZD } from './modules/vcftools_bed.nf'
include { vcftools_bed as vcftools_bedFD } from './modules/vcftools_bed.nf'
include { vcftools_bed as vcftools_bedxG } from './modules/vcftools_bed.nf'
include { vcftools_bed as vcftools_bedNC } from './modules/vcftools_bed.nf'
include { PopGenome as PopGenomeZD } from './modules/PopGenome.nf'
include { PopGenome as PopGenomeFD } from './modules/PopGenome.nf'
include { PopGenome as PopGenomexG } from './modules/PopGenome.nf'
include { PopGenome as PopGenomeAG } from './modules/PopGenome.nf'
include { PopGenome as PopGenomeNC } from './modules/PopGenome.nf'

workflow {
	contig_ch=Channel.fromPath(params.metadata)
		.splitCsv()
		.map {row -> tuple(row[0], row[1], row[2])}
		.unique()


        ZDed=ZD_BED(contig_ch)
	FDed=FD_BED(contig_ch)
        //Xged=Xgene_BED(contig_ch)
	//NCed=NCgene_BED(contig_ch)

	VCF_filtered=vcftools_filt(contig_ch)
	
	//VCF_nobeded=tabix(VCF_filtered)
	
	VCF_bededZ=vcftools_bedZD(ZDed.combine(VCF_filtered, by: [0,1,2]))
        VCF_bededD=vcftools_bedFD(FDed.combine(VCF_filtered, by: [0,1,2]))
        //VCF_bededxG=vcftools_bedxG(Xged.combine(VCF_filtered, by: [0,1,2]))
        //VCF_bededNC=vcftools_bedNC(NCed.combine(VCF_filtered, by: [0,1,2]))

	PopGenomedZ=PopGenomeZD(VCF_bededZ)	
        PopGenomedF=PopGenomeFD(VCF_bededD)
        //PopGenomedxG=PopGenomexG(VCF_bededxG)
	//PopGenomeAG=PopGenomeAG(VCF_nobeded)
        //PopGenomeNC=PopGenomeNC(VCF_bededNC)
 
}




