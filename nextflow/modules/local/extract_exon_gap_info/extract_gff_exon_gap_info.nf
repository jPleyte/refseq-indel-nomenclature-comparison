/*
 * Find all the cigar strings in a GFF that indicate a gap / difference between refseq transcript and reference genome 
 */
process EXTRACT_GFF_EXON_GAP_INFO {
    publishDir "${params.outdir}/gap", mode: 'symlink'

    input:
        path ncbi_refseq_gff_db // eg GCF_000001405.25_GRCh37.p13_genomic.gff.db
    
    output:
    path "gff_exon_gap_info.parquet", emit: gff_exon_gap_info_parquet
    path "gff_exon_gap_info.csv", emit: gff_exon_gap_info_csv

    script:
    """
    python -m rinc.io.gffutils_helper createAccessionIndex \
           --gff_db ${ncbi_refseq_gff_db} \
           --out_parquet gff_exon_gap_info.parquet \
           --out_csv gff_exon_gap_info.csv
    """
}