module load Nextflow/22.04.0

nextflow run /users/bop20pp/personal_git/SupergeneEvolutionZebraFinch/PopGenome/40Zmain.nf \
	--metadata /users/bop20pp/personal_git/SupergeneEvolutionZebraFinch/PopGenome/metadatass.csv \
	--degen /mnt/parscratch/users/bop20pp/SGPG/degeneracy-all-sites.bed \
	--gff /mnt/parscratch/users/bop20pp/SGPG/genome_files/intron_GCF_003957565.2_bTaeGut1.4.pri_genomic.gff \
	--vcfs /mnt/parscratch/users/bop20pp/SGPG/vcfs/ \
	--sample_info /users/bop20pp/personal_git/SupergeneEvolutionZebraFinch/PopGenome/sample_info.csv \
	-resume
