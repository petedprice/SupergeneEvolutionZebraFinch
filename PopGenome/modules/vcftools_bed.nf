process vcftools_bed {
    cpus = 8
    memory = '64 GB'
    time = '4h'

    label 'vcftools'

    input:
    tuple val(vcf), val(contig), val(ctg_len), file(bed), file("${contig}_filt.vcf.gz")

    output:
    tuple val(contig), val(ctg_len), env(bedtype), file('*type*vcf.gz'), file('*type*vcf.gz.tbi')    
    
    script:
    """
    #!/bin/bash
    
    basename=\$(basename $bed .bed)
    bedtype="\${basename##*_}"

    echo \$bedtype

    vcftools --gzvcf ${contig}_filt.vcf.gz --bed $bed --recode --keep-INFO-all --stdout | bgzip -c > ${contig}_filt_type\$bedtype.vcf.gz	
    tabix -p vcf ${contig}_filt_type\$bedtype.vcf.gz

    """
}
