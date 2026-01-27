#!/usr/bin/env nextflow

/*
 * Read variants look them up in CGD export, and wrote matches to csv
 */
process writeCgdNomenclatureToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path cgd_export_db
    path variants_csv
    
    output:
    path "cgd_nomenclature.csv", emit: cgd_nomenclature

    script:
    """
    python -m rinc.cgd.cgd_nomenclature \
           --cgd_db ${cgd_export_db} \
           --variants_input ${variants_csv} \
           --out cgd_nomenclature.csv
    """
}
