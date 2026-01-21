#!/usr/bin/env nextflow

/*
 * Generate the c. and p. nomenclature using the hgvs python package 
 */
process writeHgvsNomenclatureToCsv {
    publishDir "${params.outdir}", mode: 'symlink'

    input:    
    path fasta
    path variants
    
    output:
    path "hgvs_nomenclature.csv", emit: hgvs_nomenclature

    script:
    """
    python -m rinc.hgvs_uta.hgvs_nomenclature \
           --fasta ${fasta} \
           --variants ${variants} \
           --out hgvs_nomenclature.csv
    """
}
