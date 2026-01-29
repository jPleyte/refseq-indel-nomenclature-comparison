/*
 * Find all the cigar strings in a UTA that indicate a gap / difference between refseq transcript and reference genome 
 */
process EXTRACT_UTA_EXON_GAP_INFO {
    publishDir "${params.outdir}/gap", mode: 'symlink'
    
    input:
    val uta_schema

    output:
    path "uta_exon_gap_info.csv", emit: uta_exon_gap_info_csv

    script:
    """
    python ${projectDir}/bin/extract_uta_exon_gap_info.py --uta_schema ${uta_schema} --out_csv uta_exon_gap_info.csv
    """
}