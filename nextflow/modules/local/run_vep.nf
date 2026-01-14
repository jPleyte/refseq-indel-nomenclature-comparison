#!/usr/bin/env nextflow

/*
 * Generate the c. and p. nomenclature using vep
 * --use_transcript_ref (The Default): VEP will use the sequence from the transcript record. 
 *                                     This is usually what you want for HGVS nomenclature because it ensures the "c." and "p." coordinates make sense.
 * --use_given_ref: This would force VEP to use your FASTA sequence instead.
 */
process runVep {
    publishDir "${params.outdir}", mode: 'symlink'

    input:    
    path variants_vcf
    path vep_fasta
    val label
    val flag

    output:
    path "vep_${label}_output.tsv", emit: vep_output
    

    script:
    assert flag == "--use_transcript_ref" || flag == "--use_given_ref" : "Invalid VEP flag: $flag. Must be --use_transcript_ref or --use_given_ref"

    """
    # Run once with --use_transcript_ref: VEP will use the sequence from the transcript record. 
    ${params.vep}                  \
    ${flag}                        \
    --offline                      \
    --merged                       \
    --dir_cache /Users/pleyte/.vep \
    --species homo_sapiens         \
    --assembly GRCh37              \
    --fasta ${params.vep_fasta}    \
    --hgvs                         \
    --hgvsg                        \
    --symbol                       \
    --numbers                      \
    --biotype                      \
    --tab                          \
    --shift_3prime 1               \
    -i ${variants_vcf}             \
    -o vep_${label}_output.tsv
    """
}
