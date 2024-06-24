process tabix {
    cpus = 8
    memory = '64 GB'
    time = '4h'

    label 'vcftools'

    tag {'tabix_' + contig}

    input:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_filt.vcf.gz")

    output:
    tuple val(contig), val(ctg_len), env(bedtype), file('*type*vcf.gz'), file('*type*vcf.gz.tbi')
    
    script:
    """
    #!/bin/bash
    bedtype=all_SNPs    

    cp ${contig}_filt.vcf.gz ${contig}_filt_type\$bedtype.vcf.gz

    gunzip -f ${contig}_filt_type\$bedtype.vcf.gz > ${contig}_filt_type\$bedtype.vcf
    bgzip -c ${contig}_filt_type\$bedtype.vcf > ${contig}_filt_type\$bedtype.vcf.gz

    tabix -p vcf ${contig}_filt_type\$bedtype.vcf.gz

    """
}
