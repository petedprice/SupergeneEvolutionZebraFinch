# How to run this pipeline for identification of allele-specific expression 


## 1. Run nextflow/genotype_main.nf to create raw unfiltered SNPs


### Pipeline 
1. Trim reads using trimmomatic
2. Index reference genome with hisat2
3.  Allign using hisat2
4.  GATK snp calling: add readgroups, mark duplicates, split N cigar reads, haplotypecaller
6. Salmon indexing and quantification quantification 

### Commands 
~/path/to/nextflow run ./nextflow/genotyping_main.nf \
	--reads reads \
	--adapter adapters.fa \
	--ref GCF_003957565.2_bTaeGut1.4.pri_genomic.fna \
	--ref_index GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.fai \
	--cdna GCF_003957565.2_bTaeGut1.4.pri_genomic.cdna \
	--metadata full_metadata.csv \
	--contigs csvs/Z.csv \
	-resume \
	-bg


## 2. Filter SNPs

1. Genotype gVCFs from (1) using
   gatk --java-options "-Xmx4g" GenotypeGVCFs \ -R ref.fasta \ -V input.g.vcf.gz \ -O output.vcf.gz
2. Pre WASP filtering of variants (script: ./var_filt/pre_wasp_Variant_filter.sh)
   Filter SNPs that within clusters of SNPs to reduce mapping bias
   Select biallelic SNPs
3. Create text-based SNP files for each sample
   WASP/mapping/extract_vcf_snps.sh 


## 3. WASP nextflow/pipeline wasp_main.nf
### Pipeline 
Very similar to pipeline in 1 with a few changes 
This section of the pipeline removes mapping bias that may emerge in analyses of allele-specific-expression
1. Trim reads using trimmomatic
2. Index reference genome with hisat2. 
3. Alignment step 1. Standard alignment with hisat2.
4. Alignment step 2.

	a. using find_intersecting_snps.py from WASP to find reads with mapping bias
    outputs: bam file of reads that don't overlap SNPs, bamfile with original reads that overlap SNPs and need to be remapped, fastq of reads with alleles flipped to remap (remove mapping bias) 

	b. remap flipped reads identified by wasp using hisat2

	c. using filter_remapped_reads.py, filter reads where remapped read maps to new location

	d. merge non-overalpping bam from (a) and filter bam from (c)

	e. Remove duplicate reads using WASP/mapping/rmdup_pe.py
6. GATK sno calling pipeline


## 4. Filter variants from wasp pipeline
script: ./var_file/Variant_filter.sh
Filter remapped reads for clusters of 5 SNPs in a window of 95 using GATK
Filter --minGQ 30 --minDP 10 using vcftools 


##### DELETE ######
Before this, you will need to have filtered you genotyped vcfs from the first run, and filtered for depth, quality (using ../var_filt/asp_Variant_filter.sh).
You will then filter the genotyped vcfs from this run using (../var_file/Variant_filter.sh) before running phaser using the second set of filtered vcfs (../phaser)
