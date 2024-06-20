process gene_BED {
    cpus = 4
    memory = '4 GB'
    time = '32h'

    label 'bedops'

    input:
    tuple val(vcf), val(contig), val(ctg_len)

    output:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_Xgene.bed")    
    
    script:
    """
    #!/bin/bash
    cat ${params.gff} | grep $contig | 
	awk -F'\t' '\$3 == "gene"' | gff2bed > ${contig}_Xgene.bed
    """
}
