#!/bin/bash

#$ -l h_rt=1:0:0

#$ -l rmem=8G

#$ -pe smp 2

#$ -wd /fastdata/bop20pp/supergene_AB_maleref/genotype/wdir

source /usr/local/extras/Genomics/.bashrc

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 \
	-jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar \
	IndexFeatureFile \
	-F $1
