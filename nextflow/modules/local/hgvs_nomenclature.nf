#!/usr/bin/env nextflow

/*
 * Generate the c. and p. nomenclature using the hgvs python package 
 */
process hgvsNomenclature {
    publishDir "${params.outdir}", mode: 'symlink'

    input:    
    path variants
    
    output:
    path "hgvs_nomenclature.csv", emit: hgvs_nomenclature

    script:
    """
    python -m rinc.hgvs_nomenclature \
           --variants ${variants} \
           --out hgvs_nomenclature.csv
    """
}
