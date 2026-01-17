#!/usr/bin/env nextflow

/*
 * Convert tfx to a csv variant list
 */
process getTfxVariants {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path tfx_json
    
    output:
    path "variants.csv", emit: variants

    script:
    """
    python -m rinc.tfx.tfx_to_variants_csv \
           --in ${tfx_json} \
           --out variants.csv \
           --no_ccds
    """

}
