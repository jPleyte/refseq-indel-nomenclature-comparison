#!/usr/bin/env nextflow

/*
 * Convert csv variant list to annovar avinput format
 */
process writeAnnovarNomenclatureToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path multianno
    
    output:
    path "annovar_nomenclature.csv", emit: annovar_nomenclature

    script:
    """
    python -m rinc.annovar.parse_annovar_multianno \
           --annovar_multianno ${multianno} \
           --out annovar_nomenclature.csv
    """
}
