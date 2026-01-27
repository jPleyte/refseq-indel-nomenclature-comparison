#!/usr/bin/env nextflow

/*
 * Extract features from VEP's tsv and write to csv 
 */
process writeVepNomenclatureToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:    
    path vep_output_tsv
    val label
    
    output:
    path "vep_${label}_nomenclature.csv", emit: vep_nomenclature

    script:
    """
    python -m rinc.vep.vep_nomenclature    \
           --label ${label}                \
           --vep_results ${vep_output_tsv}  \
           --out vep_${label}_nomenclature.csv
    """
}
