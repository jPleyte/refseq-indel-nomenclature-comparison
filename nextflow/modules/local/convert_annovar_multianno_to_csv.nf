#!/usr/bin/env nextflow

/*
 * Convert csv variant list to annovar avinput format
 */
process convertAnnovarMultiannoToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path multianno
    
    output:
    path "annovar.avinput", emit: annovar_nomenclature

    script:
    """
    python -m rinc.etl.annovar_multianno \
           --annovar_multianno ${multianno} \
           --out annovar_nomenclature.csv
    """
}
