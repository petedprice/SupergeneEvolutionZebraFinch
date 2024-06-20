process PopGenome {
    cpus = 16
    memory = '64 GB'
    time = '4h'
    label 'R'

    publishDir 'PopGenome_Stats', mode: 'copy', overwrite: true, pattern: '*_PopGenome_Stats.csv'

    input:
    tuple val(contig), val(ctg_len), val(bedtype), file(vcf), file(index)

    output:
    tuple val(contig), val(bedtype), file("${bedtype}_${contig}_PopGenome_Stats.csv")

    script:
    """
    #!/bin/bash

    Rscript ${projectDir}/scripts/PopGenome.R $vcf $contig $ctg_len ${params.sample_info} $bedtype

    
    """
}
