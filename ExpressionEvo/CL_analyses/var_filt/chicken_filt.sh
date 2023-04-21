#!/bin/bash

#$ -l h_rt=1:0:0

#$ -l rmem=8G

#$ -pe smp 2

#$ -wd /fastdata/bop20pp/supergene_AB_maleref/chicken/genotyped_vcfs/wdir

source /usr/local/extras/Genomics/.bashrc

#$ -P ressexcon
#$ -q ressexcon.q

out=${1%".genotyped.vcf.gz"}
file=${out#"/fastdata/bop20pp/supergene_AB_maleref/chicken/genotyped_vcfs/"}

echo $out

tabix -p vcf $1

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar VariantFiltration \
	-R /fastdata/bop20pp/supergene_AB_maleref/chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa \
	-V $1 \
	-O ${out}.tempclus.vcf \
	-cluster 5 \
	-window 100


#cp $1 $out.tempclus.filt.vcf.gz
#gunzip $out.tempclus.filt.vcf.gz
bcftools view -f PASS $out.tempclus.vcf > $out.tempclus.filt.vcf


gunzip $out.tempclus.filt.vcf.gz

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 \
	-jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar \
	SelectVariants \
	-R /fastdata/bop20pp/supergene_AB_maleref/chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa \
	-V $out.tempclus.filt.vcf \
	--restrict-alleles-to BIALLELIC \
	--select-type-to-include SNP \
	-O $out.tempclus.filt.ba.vcf


vcftools --vcf $out.tempclus.filt.ba.vcf --out /fastdata/bop20pp/supergene_AB_maleref/chicken/filtered_genotyped_vcfs/${file}_filtered --minGQ 30 --minDP 10 --recode --recode-INFO-all

