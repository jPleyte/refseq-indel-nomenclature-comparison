#!/usr/bin/env nextflow

/*
 * Read variants look them up in CGD export, and wrote matches to csv
 */
process writeCgdNomenclatureToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path variants_csv
    
    output:
    path "cgd_nomenclature.csv", emit: cgd_nomenclature

    script:
    """
    python -m rinc.cgd.cgd_nomenclature \
           --cgd_db ${params.cgd_export_csv} \
           --variants_input ${variants_csv} \
           --out cgd_nomenclature.csv
    """
}
