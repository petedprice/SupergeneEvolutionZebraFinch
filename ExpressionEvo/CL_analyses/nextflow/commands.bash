~/software/nextflow run /home/bop20pp/software/supergene_wogtf/nextflow/genotyping_main.nf \
	--reads reads \
	--adapter adapters.fa \
	--ref GCF_003957565.2_bTaeGut1.4.pri_genomic.fna \
	--ref_index GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.fai \
	--cdna GCF_003957565.2_bTaeGut1.4.pri_genomic.cdna \
	--metadata full_metadata.csv \
	--contigs csvs/Z.csv \
	-resume \
	-bg
