#!/usr/bin/env nextflow

/*
 * Convert csv variant list to annovar avinput format
 */
process csvToAvinput {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path variant_csv_file
    
    output:
    path "annovar.avinput", emit: annovar_avinput

    script:
    """
    python -m rinc.etl.csv_to_avinput \
           --in ${variant_csv_file} \
           --out annovar.avinput
    """
}
