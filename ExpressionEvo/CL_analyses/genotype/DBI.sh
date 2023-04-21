#!/bin/bash

#$ -l h_rt=2:0:0

#$ -l rmem=4G

#$ -pe smp 1

#$ -wd /fastdata/bop20pp/supergene_AB_maleref/genotype/wdir

source /usr/local/extras/Genomics/.bashrc

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 \
	-jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar \
	GenomicsDBImport \
	-V $1 \
	-V $2 \
	-V $3 \
	-V $4 \
	-V $5 \
	-V $6 \
	-V $7 \
	-V $8 \
	-V $9 \
	-V ${10} \
	--genomicsdb-workspace-path ../Z_database \
	--intervals ${11} \
	-R /fastdata/bop20pp/supergene_AB_maleref/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna
