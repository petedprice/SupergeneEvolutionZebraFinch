process vcftools_filt {
    cpus = 1
    memory = '4 GB'
    time = '4h'

    label 'vcftools'

    tag {'vcftools_filt_' + contig }

    input:
    tuple val(vcf), val(contig), val(ctg_len)

    output:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_filt.vcf.gz")    
    
    script:
    """
    #!/bin/bash

    vcftools --gzvcf ${params.vcfs}/$vcf --recode --minDP 5 --minGQ 20 --max-missing 0.5 --stdout | gzip -c > ${contig}_filt.vcf.gz	

    """
}
