/*
 * Merge the UTA and GFF lists of transcritps that have known gaps with reference genome
 */
process COMBINE_GFF_AND_UTA_EXON_GAP_INFO {
    publishDir "${params.outdir}/gap", mode: 'symlink'
    
    input:
    path gff_exon_gap_info_parquet
    path uta_exon_gap_info_csv

    output:
    path "gff_and_uta_exon_gap_info.csv", emit: gff_and_uta_exon_gap_info

    script:
    """
    python -m rinc.analysis.combine_gff_and_uta_gap_info \
           --gff_gaps ${gff_exon_gap_info_parquet}          \
           --uta_gaps ${uta_exon_gap_info_csv}              \
           --out_csv gff_and_uta_exon_gap_info.csv
    """
}