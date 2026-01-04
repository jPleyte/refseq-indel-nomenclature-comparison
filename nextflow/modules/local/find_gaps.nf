#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    output:
        path "gap_transcript_variants.csv", emit: gap_transcript_variants

    script:
    """
    python3 ${projectDir}/python/src/rinc/etc/find_gaps.py --out gap_transcript_variants.csv --fasta ${projectDir}/Homo_sapiens_assembly19.fasta
    """
}
