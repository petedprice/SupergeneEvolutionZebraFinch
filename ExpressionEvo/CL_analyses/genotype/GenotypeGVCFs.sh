#!/bin/bash

#$ -l h_rt=2:0:0

#$ -l rmem=4G

#$ -pe smp 1

#$ -wd /fastdata/bop20pp/supergene_AB_maleref/genotype/wdir

source /usr/local/extras/Genomics/.bashrc

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 \
	-jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar \
	GenotypeGVCFs \
	-R $1 \
	-V gendb://$2 \
	-O /fastdata/bop20pp/supergene_AB_maleref/genotype/Z.vcf

