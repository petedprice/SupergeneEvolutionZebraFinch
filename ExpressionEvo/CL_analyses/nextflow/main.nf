reads_ch = Channel.fromFilePairs(params.reads + '/*{1,2}.fastq.gz')
ref=file(params.ref)
ref_index=file(params.ref_index)
metadata=file(params.metadata)
adapter=file(params.adapter)
contigs=file(params.contigs)

//targetsnps=params.targetsnps

ref_ch=channel
    .fromPath(metadata)
    .splitCsv()
    .map {row ->tuple(row[0],row[1])}


contig_ch=Channel
    .from(contigs)
    .splitCsv()
    .view()
  


process trim {
    errorStrategy 'finish'

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = 4
    memory = '16 GB'
    time = '4h'
    
    tag {'trim_' + species + '_' + sid }


    publishDir 'paired', mode: 'copy', overwrite: true, pattern: '*paired*'

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
    cache 'lenient'

    errorStrategy 'finish'
    queue = "ressexcon.q"
    cpus = 16
    memory = '192 GB'
    time = '1h'
    clusterOptions = { '-P ressexcon' }

    
    tag {'hisat index'}



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
    cache 'lenient'
    errorStrategy 'finish'    
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
    cache 'lenient' 
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
    cache 'lenient' 
    publishDir 'bam', mode: 'copy', overwrite: true, pattern: '*bam'
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
    cache 'lenient'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    tag {'SplitNCigarReads_' + '_' + sid + '_' + contig }

    publishDir 'bam', mode: 'copy', overwrite: true, pattern: '*bam*'

    queue = "molecosh.q"
    clusterOptions = { '-P molecosh' }

    cpus = { 4 * task.attempt }
    memory = { 12.GB * task.attempt }
    time = { 12.hour }
    maxRetries 10


    input:
    tuple val(sid), file("${sid}_md.RG.bam"), val(contig) from  mark_dup1.combine(contig_ch)
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index

    output:
    tuple val(sid), file("${contig}_${sid}_md_cig.RG.bam"), file("${contig}_${sid}_md_cig.RG.bam.bai"), val(contig), file('ref.dict') into cigar1, cigar2

    """
    #!/bin/bash
    module load apps/python/anaconda3-4.2.0
    source activate gatk
    mkdir tmp
    gatk CreateSequenceDictionary -R ref.fasta -O ref.dict 
    samtools index ${sid}_md.RG.bam
    gatk SplitNCigarReads \
	--java-options "-Xmx${task.memory.toMega()}M -XX:ParallelGCThreads=${task.cpus}" \
	-R ref.fasta -I ${sid}_md.RG.bam \
	-O ${contig}_${sid}_md_cig.RG.bam \
	-L $contig \
	--tmp-dir tmp 

    samtools index ${contig}_${sid}_md_cig.RG.bam
    """

}


process snp_calling_gvcf {
    cache 'lenient'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }

    tag {'snp_calling_gvcf_' + '_' + sid + '_' + contig}
    publishDir 'gvcfs', mode: 'copy', overwrite: true, pattern: '*.g.vcf.gz'

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = { 2 * task.attempt }
    memory = { 192.GB }
    time = { 12.hour }
    maxRetries 10


    input:
    tuple val(sid), file("${contig}_${sid}_md_cig.RG.bam"), file("${contig}_${sid}_md_cig.RG.bam.bai"), val(contig), file('ref.dict') from cigar2
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index

    output:
    tuple val(sid), file("${contig}_${sid}_md_cig.RG.bam"), file("${contig}_${sid}.g.vcf.gz"), val(contig) into snp_gvcf_called1
    tuple file("${contig}_${sid}.g.vcf.gz"), val(sid), val(contig)  into to_genotype


    script:
    """
    #!/bin/bash
    module load apps/python/anaconda3-4.2.0
    source activate gatk
    mkdir tmp
    gatk HaplotypeCaller \
        --java-options "-Xmx${task.memory.toMega()}M -XX:ParallelGCThreads=${task.cpus}" \
        -R ref.fasta \
        -I ${contig}_${sid}_md_cig.RG.bam \
        -O ${contig}_${sid}.g.vcf.gz \
        -native-pair-hmm-threads ${task.cpus} \
	-L $contig \
        -ERC GVCF \
       	--tmp-dir tmp
    rm -r tmp
    """

}

/*
process genotpye_gvcfs {
    errorStrategy 'finish'

    tag {'genotpye_gvcfs'}
    publishDir 'genotyped_vcfs', mode: 'copy', overwrite: true, pattern: '*vcf.gz'

    queue = "ressexcon.q"
    clusterOptions = { '-P ressexcon' }
    cpus = 8
    memory = '64 GB'
    time = '12h'

    input:
    tuple file("${contig}_${sid}.g.vcf.gz"), val(sid) from to_genotype.collect()
    file('ref.fasta') from ref
    file('ref.fasta.fai') from ref_index

    output:
    file("genotyped.vcf") into genotyped1

    script:
    """
    #!/bin/bash
    
    source /usr/local/extras/Genomics/.bashrc
    module load apps/gatk/4.1.4/binary
    mkdir tmp

    bcftools concat *.vcf.gz > ${sid}.g.vcf
    gzip ${sid}.g.vcf    
    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /usr/local/community/Genomics/apps/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar GenotypeGVCFs \
	-R ref.fasta \
	-V ${sid}.g.vcf.gz \
	-O ${sid}.genotyped.vcf.gz
    
    """

}

*/
