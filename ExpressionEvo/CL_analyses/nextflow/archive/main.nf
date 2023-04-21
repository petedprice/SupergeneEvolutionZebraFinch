reads_ch = Channel.fromFilePairs(params.reads + '/*{1,2}.fastq.gz')
ref=file(params.ref)
ref_index=file(params.ref_index)
metadata=file(params.metadata)
adapter=file(params.adapter)
//targetsnps=params.targetsnps

ref_ch=channel
    .fromPath(metadata)
    .splitCsv()
    .map {row ->tuple(row[0],row[1])}


process trim {


    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = 4
    memory = '8 GB'
    time = '2h'
    
    tag {'trim_' + species + '_' + sid }


    //publishDir 'paired', mode: 'copy', overwrite: true, pattern: '*paired*'

    input:
    tuple val(sid), file(reads), val(species) from reads_ch.combine(ref_ch, by:0)

    output:
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") into trimmed1
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") into trimmed2

    script:

    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    trimmomatic PE -phred33 $reads ${species}_${sid}_forward_paired.fastq.gz ${sid}_forward_unpaired.fastq.gz ${species}_${sid}_reverse_paired.fastq.gz ${sid}_reverse_unpaired.fastq.gz ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:95
    """
}



process index_hisat2 {

    queue = "ressexcon.q"
    cpus = 16
    memory = '192 GB'
    time = '1h'
    clusterOptions = { '-P ressexcon' }

    
    tag {'hisat index'}


    publishDir 'index', mode: 'copy', overwrite: true, pattern: '*'

    output:
    file('hisat_index') into hisat_indexed


    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    hisat2-build $ref index -p 16
    mkdir hisat_index
    mv index*ht2 hisat_index
    """



}


process allignment_hisat2 {
    
    tag {'allign_' + '_' + sid }

    queue = "ressexcon.q"
    cpus = 16
    memory = '8 GB'
    time = '4h'
    clusterOptions = { '-P ressexcon' }
    

    //publishDir 'bam', mode: 'copy', overwrite: true, pattern: '*bam'

    input:
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") from trimmed2
    file('hisat_index') from hisat_indexed

    output:
    tuple val(sid), file("${sid}.bam") into alligned1
    tuple val(sid), file("${sid}.bam") into alligned2     

    script: 
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    hisat2 \
	-x hisat_index/index \
	-1 ${species}_${sid}_forward_paired.fastq.gz \
	-2 ${species}_${sid}_reverse_paired.fastq.gz \
        --summary-file ${sid}_hisat2.summary.log \
        --threads 16 \
            | samtools view -bS -F 4 -F 8 -F 256 - > ${sid}.bam
    """

}

process RG_add {
    errorStrategy 'finish'

    tag {'RG_' + '_' + sid }

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }

    cpus = 4
    memory = '8 GB'
    time = '3h'

    input:
    tuple val(sid), file("${sid}.bam") from alligned1

    output:
    tuple val(sid), file("${sid}.RG.bam") into readgrouped1
    tuple val(sid), file("${sid}.RG.bam") into readgrouped2

    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    source activate gatk
    mkdir `pwd`/tmp
    picard AddOrReplaceReadGroups I=${sid}.bam O=${sid}.RG.bam SORT_ORDER=coordinate RGID=zf RGLB=zf RGPL=illumina RGSM=zf RGPU=zf CREATE_INDEX=True TMP_DIR=`pwd`/tmp
    """

}

process mark_duplicates {
    errorStrategy 'finish'

    tag {'mark_duplicates_' + '_' + sid }

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }

    cpus = 4
    memory = '8 GB'
    time = '2h'

    input:
    tuple val(sid), file("${sid}.RG.bam") from readgrouped1

    output:
    tuple val(sid), file("${sid}_md.RG.bam") into mark_dup1

    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    mkdir tmp
    samtools sort ${sid}.RG.bam -o ${sid}.RG.sorted.bam
    picard MarkDuplicates \
      	I=${sid}.RG.sorted.bam \
      	O=${sid}_md.RG.bam \
      	M=${sid}_marked_dup_metrics.txt \
	TMP_DIR=tmp
    """

}

process SplitNCigarReads {
    errorStrategy 'finish'

    tag {'SplitNCigarReads_' + '_' + sid }

    publishDir 'bam', mode: 'copy', overwrite: true, pattern: '*bam'

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }

    cpus = 8
    memory = '32 GB'
    time = '2h'

    input:
    tuple val(sid), file("${sid}_md.RG.bam") from  mark_dup1
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index
    
    output:
    tuple val(sid), file("${sid}_md_cig.RG.bam") into cigar1, cigar2

    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    module load apps/gatk/4.1.4/binary
    mkdir tmp
    picard CreateSequenceDictionary R=ref.fasta O=ref.dict 
    samtools index ${sid}_md.RG.bam
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx32g -jar /usr/local/packages/apps/gatk/4.1.4/binary/gatk-package-4.1.4.0-local.jar SplitNCigarReads -R ref.fasta -I ${sid}_md.RG.bam -O ${sid}_md_cig.RG.bam -L CM020861.1 --tmp-dir tmp 
    """

}


/*
process snp_calling {
    errorStrategy 'finish'

    tag {'haplotypecaller_' + '_' + sid }
    publishDir 'vcfs', mode: 'copy', overwrite: true, pattern: '*vcf.gz'

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = 16
    memory = '192 GB'
    time = '24h'

    input:
    tuple val(sid), file("${sid}_md_cig.RG.bam") from cigar1
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index

    output:
    tuple val(sid), file("${sid}_md_cig.RG.bam"), file("${sid}.vcf.gz") into snp_called1

    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    module load apps/gatk/4.1.4/binary
    source activate picard 
    mkdir tmp
    samtools faidx ref.fasta
    samtools index ${sid}_md_cig.RG.bam
    picard CreateSequenceDictionary R=ref.fasta O=ref.dict

    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx192g -XX:ParallelGCThreads=16 -jar /usr/local/packages/apps/gatk/4.1.4/binary/gatk-package-4.1.4.0-local.jar HaplotypeCaller \
       	-R ref.fasta \
        -I ${sid}_md_cig.RG.bam \
        -O ${sid}.vcf.gz \
	-L CM020861.1 \
        -native-pair-hmm-threads 8 \
       	--tmp-dir tmp
    """

}



process snpfilter {

    errorStrategy 'finish'
    tag {'featurecount_' + '_' + sid }
    publishDir 'vcfs', mode: 'copy', overwrite: true, pattern: '*vcf.gz'


    input:
    tuple val(sid), file("${sid}_md_cig.RG.bam"), file("${sid}.vcf.gz") from snp_called1
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index

    output:
    tuple val(sid), file("${sid}_final.vcf"), file("${sid}_md_cig.RG.bam") into filtered1
    file("${sid}_marked.vcf.gz") into unfiltered1
    file("${sid}_final.vcf") into filtered2


    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    tabix -p vcf ${sid}.vcf.gz
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar VariantFiltration \
	-R ref.fasta \
	-V ${sid}.vcf.gz \
	-O clus.vcf.gz \
	-cluster 10 \
	-window 145 
    ##CLUSTERS MARKED
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar VariantFiltration \
	-V clus.vcf.gz \
	-O ${sid}_marked.vcf.gz \
	--filter-expression "DP < 10" \
	--filter-name "depth_filter"
    #DEPTH FILTERED 
    bcftools view -i 'FILTER="PASS"' depth.vcf.gz > ${sid}_final.vcf
    
    #REMOVED NON-PASS VARIANTS
    
    #FS > 30
    # QD < 2
    # DP < 10
    # GQ < 30
    
    

    
    """

}


process ASEReadCounter { 

    errorStrategy 'finish'
    tag {'featurecount_' + '_' + sid }
    
    input:
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index
    tuple val(sid), file("${sid}_final.vcf"), file("${sid}_md_cig.RG.bam") from filtered1

    output:
    

    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    source activate gatk
    gatk --java-options "-Xmx4g" -R ref.fasta -I ${sid}_md_cig.RG.bam -V ${sid}_final.vcf -O ASE_${sid}.txt

    """

}

*/

process snp_calling_gvcf {
    errorStrategy 'finish'

    tag {'snp_calling_gvcf_' + '_' + sid }
    publishDir 'gvcfs', mode: 'copy', overwrite: true, pattern: '*.g.vcf.gz'

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = 16
    memory = '32 GB'
    time = '12h'

    input:
    tuple val(sid), file("${sid}_md_cig.RG.bam") from cigar2
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index

    output:
    tuple val(sid), file("${sid}_md_cig.RG.bam"), file("${sid}.g.vcf.gz") into snp_gvcf_called1
    file("${sid}.g.vcf.gz") into to_genotype


    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    module load apps/gatk/4.1.4/binary
    source activate picard
    mkdir tmp
    samtools faidx ref.fasta
    samtools index ${sid}_md_cig.RG.bam
    picard CreateSequenceDictionary R=ref.fasta O=ref.dict

    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx192g -XX:ParallelGCThreads=16 -jar /usr/local/packages/apps/gatk/4.1.4/binary/gatk-package-4.1.4.0-local.jar HaplotypeCaller \
        -R ref.fasta \
        -I ${sid}_md_cig.RG.bam \
        -O ${sid}.g.vcf.gz \
        -L CM020861.1 \
        -native-pair-hmm-threads 8 \
        -ERC GVCF \
       	--tmp-dir tmp
    rm -r tmp
    #echo ${sid} > sample.txt
    #bcftools reheader -s sample.txt ${sid}.g.vcf.gz -o ${sid}.g.vcf.gz
    """

}


process genotpye_gvcfs {
    errorStrategy 'finish'

    tag {'genotpye_gvcfs'}
    publishDir 'genotyped_vcfs', mode: 'copy', overwrite: true, pattern: '*vcf.gz'

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = 4
    memory = '8 GB'
    time = '2h'

    input:
    file(samples) from to_genotype.collect()
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index

    output:
    file("genotyped.vcf") into genotyped1

    script:
    """
    #!/bin/bash
    
    source /usr/local/extras/Genomics/.bashrc
    module load apps/gatk/4.1.4/binary
    

    for i in *.g.vcf.gz
    do
    echo $i > sample.txt 
    bcftools reheader -s sample.txt $i -o rh_$i
    done


    mkdir tmp
    echo -n java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar GenomicsDBImport \\ >> command.bash 
    for i in rh*.g.vcf.gz
    do
    echo -n -V \$i \\ >> command.bash
    done
    echo --genomicsdb-workspace-path my_database --tmp-dir=tmp -L CM020861.1 >> command.bash
    bash command.bash
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar GenotypeGVCFs \
	-R ref.fasta \
	-V gendb://my_database \
	-O genotyped.vcf
    
    """

}

