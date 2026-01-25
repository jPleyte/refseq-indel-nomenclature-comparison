#!/usr/bin/env nextflow

/*
 * Write exon coordinates and cigar strings for each transcripts 
 */
process writeExonDetail {
    input:
    path ncbi_refseq_gff_db
    path ncbi_refseq_gff_accession_index_df
    val tool_tuples

    output:
    path "transcript_exon_detail.csv", emit: transcript_exon_detail

    script:
    // Iterate over the list to build a string like: --tool annovar file1.csv --tool vep file2.csv
    def tool_args = tool_tuples.collect { _label, file -> "--transcripts ${file}" }.join(' ')
    """
    python -m rinc.util.write_exon_detail ${tool_args} \
              --gff_db ${ncbi_refseq_gff_db} \
              --accession_index ${ncbi_refseq_gff_accession_index_df} \
              --out transcript_exon_detail.csv
    """
}