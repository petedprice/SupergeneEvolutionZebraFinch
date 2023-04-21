
#!/bin/bash

#$ -l h_rt=2:0:0

#$ -l rmem=8G

#$ -pe smp 4

#$ -wd /fastdata/bop20pp/supergene_AB_maleref/phaser/wdir

python ~/software/phaser/phaser/phaser.py \
	--vcf  /fastdata/bop20pp/supergene_AB_maleref/filtered_genotyped_vcfs/${1}_tgu_Sample_${2}_filtered.recode.vcf.gz \
	--bam /fastdata/bop20pp/supergene_AB_maleref/bam/${1}_tgu_Sample_${2}_md_cig.RG.bam \
	--paired_end 1 --mapq 60 --baseq 10 \
	--sample zf --threads 4 --o /fastdata/bop20pp/supergene_AB_maleref/phaser/$2/${1}_${2}_phaser \
	--id_separator +
