#!/usr/bin/env nextflow

/*
 * Read Tfx json and write out a nomenclature csv.
 */
process writeTfxNomenclatureToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path fasta
    path tfx_variants_file

    output:
    path "tfx_nomenclature.csv", emit: tfx_nomenclature

    script:
    """
    python -m rinc.tfx.tfx_nomenclature \
        --fasta ${fasta} \
        --tfx_input ${tfx_variants_file} \
        --out tfx_nomenclature.csv
    """
}    
