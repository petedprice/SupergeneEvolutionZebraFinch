#!/bin/bash

#$ -l h_rt=2:0:0

#$ -l rmem=16G

#$ -pe smp 4

#$ -wd /fastdata/bop20pp/supergene_AB_maleref/phaser/wdir

source /usr/local/extras/Genomics/.bashrc
module load apps/python/anaconda2-4.2.0
source activate phaser

mkdir /fastdata/bop20pp/supergene_AB_maleref/phaser/$2/ASE

python2.7 ~/software/phaser/phaser_gene_ae/phaser_gene_ae.py \
	--haplotypic_counts /fastdata/bop20pp/supergene_AB_maleref/phaser/$2/${1}_${2}_phaser.haplotypic_counts.txt \
	--features  /fastdata/bop20pp/supergene_AB_maleref/phaser/beds/${1}.bed \
	--o /fastdata/bop20pp/supergene_AB_maleref/phaser/$2/ASE/${1}_${2}.phaserASE.txt \
	--id_separator +
