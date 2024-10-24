nextflow.enable.dsl=2

include { FD_BED } from './modules/FD_BED.nf'
include { ZD_BED } from './modules/ZD_BED.nf'
include { Xgene_BED } from './modules/Xgene_BED.nf'
include { NCgene_BED } from './modules/NCgene_BED.nf'
include { NCXG_BED } from './modules/NCXG_BED.nf'

include { tabix } from './modules/tabix.nf'
include { vcftools_filt } from './modules/vcftools_filt.nf'

include { pseudohaploid } from './modules/pseudohaploid.nf'

include { vcftools_bed as vcftools_bedZD } from './modules/vcftools_bed.nf'
include { vcftools_bed as vcftools_bedFD } from './modules/vcftools_bed.nf'
include { vcftools_bed as vcftools_bedxG } from './modules/vcftools_bed.nf'
include { vcftools_bed as vcftools_bedNC } from './modules/vcftools_bed.nf'
include { vcftools_bed as vcftools_bedNCXG } from './modules/vcftools_bed.nf'

include { PH_PopGenome as PopGenomeZD } from './modules/PH_PopGenome.nf'
include { PH_PopGenome as PopGenomeFD } from './modules/PH_PopGenome.nf'
include { PH_PopGenome as PopGenomexG } from './modules/PH_PopGenome.nf'
include { PH_PopGenome as PopGenomeAG } from './modules/PH_PopGenome.nf'
include { PH_PopGenome as PopGenomeNC } from './modules/PH_PopGenome.nf'
include { PH_PopGenome as PopGenomeNCXG } from './modules/PH_PopGenome.nf'

workflow {
	contig_ch=Channel.fromPath(params.metadata)
		.splitCsv()
		.map {row -> tuple(row[0], row[1], row[2])}
		.unique()


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
}




