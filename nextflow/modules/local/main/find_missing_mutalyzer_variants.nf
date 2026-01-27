#!/usr/bin/env nextflow

/*
 * Read the variant transcripts file, identify the variants that aren't in our local Mutalyzer cache
 * and write out batch files for sumisstion to Mutalyzer Batch Processor.
 */
process find_missing_mutalyzer_variants {
    publishDir "${params.outdir}/mutalyzer_batch/", mode: 'symlink'

    input:
    path variant_transcripts
    path mutalyzer_cache

    output:
    path "*.txt", emit: batch_files, optional: true

    script:
    """
    python -m rinc.mutalyzer.final_csv_to_mutalyzer_batch \
        --variant_nomenclature ${variant_transcripts} \
        --mutalyzer_cache ${mutalyzer_cache} \
        --out_dir .
    """
}