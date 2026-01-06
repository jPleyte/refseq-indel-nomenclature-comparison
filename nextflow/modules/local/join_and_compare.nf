#!/usr/bin/env nextflow

/*
 * Gather gap, hgvs, and annovar information; perform comparison, and write results to csv.
 */
process joinAndCompare {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path gaps_and_variants
    path hgvs_nomenclature    
    path annovar_nomenclature    
    
    output:
    path "gap_nomenclature_comparison.csv", emit: gap_nomenclature_comparison

    script:
    """
    python -m rinc.join_and_compare \
        --gap_variants ${gaps_and_variants} \
        --hgvs_nomenclature ${hgvs_nomenclature} \
        --annovar_nomenclature ${annovar_nomenclature} \
        --out gap_nomenclature_comparison.csv
    """
}
