#!/usr/bin/env nextflow

/*
 * Extract features from SnpEff's tsv and write to csv 
 */
process convertSnpEffToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path snpeff_tsv
    
    output:
    path "snpeff_nomenclature.csv", emit: snpeff_nomenclature

    script:
    """
    python -m rinc.etl.process_snpeff \
           --in ${snpeff_tsv} \
           --out snpeff_nomenclature.csv
    """
}
