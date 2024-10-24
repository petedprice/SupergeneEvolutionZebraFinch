module load Nextflow/22.04.0

nextflow run /users/bop20pp/personal_git/SupergeneEvolutionZebraFinch/PopGenome/PH_main.nf \
	--metadata /users/bop20pp/personal_git/SupergeneEvolutionZebraFinch/PopGenome/metadata_full.csv \
	--degen /mnt/parscratch/users/bop20pp/SGPG/degeneracy-all-sites.bed \
	--vcfs /mnt/parscratch/users/bop20pp/SGPG/vcfs/ \
	--sample_info /users/bop20pp/personal_git/SupergeneEvolutionZebraFinch/PopGenome/sample_info.csv \
	-resume
