#!/usr/bin/env nextflow

/*
 * Query UTA database to find siple insertion or deletions and then create downstream variants.
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
