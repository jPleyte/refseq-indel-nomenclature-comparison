#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process findGapVariants {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    val uta_schema
    path fasta
    

    output:
    path "gaps_and_variants.csv", emit: gaps_and_variants

    script:
    """
    python -m rinc.etl.find_gap_variants \
           --uta_schema ${uta_schema} \
           --fasta ${fasta} \
           --out gaps_and_variants.csv
    """
}
