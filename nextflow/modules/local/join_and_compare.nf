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
    path snpeff_nomenclature
    path mutalyzer_nomenclature
    path vep_refseq_nomenclautre
    path vep_hg19_nomenclature

    output:
    path "gap_nomenclature_all_transcripts_all_fields.csv", emit: gap_nomenclature_all_transcripts_all_fields

    script:
    """
    python -m rinc.join_and_compare \
        --gap_variants ${gaps_and_variants} \
        --hgvs_nomenclature ${hgvs_nomenclature} \
        --annovar_nomenclature ${annovar_nomenclature} \
        --snpeff_nomenclature ${snpeff_nomenclature} \
        --mutalyzer_nomenclature ${mutalyzer_nomenclature} \
        --vep_refseq_nomenclautre ${vep_refseq_nomenclautre} \
        --vep_hg19_nomenclature ${vep_hg19_nomenclature} \
        --out gap_nomenclature_all_transcripts_all_fields.csv
    """
}
