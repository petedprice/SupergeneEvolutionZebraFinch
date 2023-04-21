#!/bin/bash

#$ -l h_rt=1:0:0

#$ -l rmem=8G

#$ -pe smp 2

#$ -wd /fastdata/bop20pp/supergene_AB_maleref/genotyped_vcfs/wdir

source /usr/local/extras/Genomics/.bashrc

out=${1%"..vcf.gz"}
file=${out#"/fastdata/bop20pp/supergene_AB_maleref/gvcfs/"}

echo $out

tabix -p vcf $1

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar VariantFiltration \
	-R /fastdata/bop20pp/supergene_AB_maleref/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna \
	-V $1 \
	-O ${out}.tempclus.vcf \
	-cluster 5 \
	-window 95


bcftools view -f PASS $out.tempclus.vcf > $out.tempclus.filt.vcf

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 \
	-jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar \
	SelectVariants \
	-R /fastdata/bop20pp/supergene_AB_maleref/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna \
	-V $out.tempclus.filt.vcf \
	--restrict-alleles-to BIALLELIC \
	--select-type-to-include SNP \
	-O $out.tempclus.filt.ba.vcf


vcftools --vcf $out.tempclus.filt.ba.vcf --out /fastdata/bop20pp/supergene_AB_maleref/filtered_genotyped_vcfs/${file}_filtered --minGQ 30 --minDP 10 --recode --recode-INFO-all

