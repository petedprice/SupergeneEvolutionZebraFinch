reads_ch = Channel.fromFilePairs(params.reads + '/*{1,2}.fastq.gz')
gtf=file(params.gtf)
ref=file(params.ref)
metadata=file(params.metadata)
adapter=file(params.adapter)
//targetsnps=params.targetsnps

ref_ch=channel
    .fromPath(metadata)
    .splitCsv()
    .map {row ->tuple(row[0],row[1])}


process trim {

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


process prep {


    output:
    tuple file('genome.ss'), file('genome.exon') into prepped

    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    hisat2_extract_splice_sites.py $gtf > genome.ss
    hisat2_extract_exons.py $gtf > genome.exon
    """

}




process index_hisat2 {

    //queue = 'ressexcon.q'
    cpus = 16
    memory = '80 GB'
    time = '4h'

    tag {'star index'}


    publishDir 'index', mode: 'copy', overwrite: true, pattern: '*'

    input:
    tuple file('genome.ss'), file('genome.exon') from prepped

    output:
    file('star_index') into star_indexed


    script:
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    #/home/bop20pp/software/STAR/source/STAR --runThreadN 6 \
#	--runMode genomeGenerate \
#	--genomeDir star_index \
#	--genomeFastaFiles $ref \
#	--sjdbGTFfile $gtf \
#	--sjdbOverhang 99 \
#	--genomeSAindexNbases 12 \
#	--limitGenomeGenerateRAM 25000000000 \
#	--sjdbGTFtagExonParentTranscript Parent
 
    hisat2-build --exon genome.exon --ss genome.ss $ref index -p 16
    """



}


process allignment_star {
    

    cpus = 4
    memory = '16 GB'
    time = '2h'

    input:
    tuple val(species), val(sid), file("${species}_${sid}_forward_paired.fastq.gz"), file("${species}_${sid}_reverse_paired.fastq.gz") from trimmed2
    file('star_index') from star_indexed

    output:
    tuple val(sid), file("star_alligned_${sid}") into alligned     

    script: 
    """
    #!/bin/bash
    source /usr/local/extras/Genomics/.bashrc
    /home/bop20pp/software/STAR/source/STAR --genomeDir star_index/ \
	--runThreadN 8 \
	--readFilesIn ${species}_${sid}_forward_paired.fastq.gz ${species}_${sid}_reverse_paired.fastq.gz \
	--outFileNamePrefix star_alligned_${sid} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard \
	--quantMode GeneCounts \
	--readFilesCommand zcat
    """

}

process snp_calling {


    input:
    tuple val(sid), file("star_alligned_${sid}") from alligned


    output:


    script:
    """
    source /usr/local/extras/Genomics/.bashrc
    gatk --java-options "-Xmx4g" HaplotypeCaller \
	-R $ref \
	-I input.bam \
	-O ${sid}.vcf.gz \
	-ERC GVCF
    """

}

/*

process allele_quant {


    input:
    tuple val(sid), file("star_alligned_${sid}") from alligned

    output:


    script:
    """
    source /usr/local/extras/Genomics/.bashrc
    gatk ASEReadCounter \ 
	-R $ref \ 
	-I ${sid}.bam \ 
	-V targetsites.vcf.gz \ 
	-O ${sid}_allele_counts.table
    """

}


*/
