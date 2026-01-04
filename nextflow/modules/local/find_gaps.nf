#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process findGaps {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fasta

    output:
    path "gap_transcript_variants.csv", emit: gap_transcript_variants

    script:
    """
    python -m rinc.etl.find_gaps --out gap_transcript_variants.csv --fasta ${fasta}
    """
}
