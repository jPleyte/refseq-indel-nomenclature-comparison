include { EXTRACT_GFF_EXON_GAP_INFO } from '../../modules/local/extract_exon_gap_info/extract_gff_exon_gap_info'
include { EXTRACT_UTA_EXON_GAP_INFO } from '../../modules/local/extract_exon_gap_info/extract_uta_exon_gap_info'
include { COMBINE_GFF_AND_UTA_EXON_GAP_INFO } from '../../modules/local/extract_exon_gap_info/combine_gff_and_uta_exon_gap_info'

workflow EXTRACT_EXON_GAP_INFO {
    take:        
        // Channel pointing to GCF_000001405.25_GRCh37.p13_genomic.gff.db
        ch_ncbi_gff_db  

        // Uta schema like uta_20240523b
        uta_schema


    main:
        EXTRACT_GFF_EXON_GAP_INFO(ch_ncbi_gff_db)
        EXTRACT_UTA_EXON_GAP_INFO(uta_schema)
        COMBINE_GFF_AND_UTA_EXON_GAP_INFO(EXTRACT_GFF_EXON_GAP_INFO.out.gff_exon_gap_info_parquet, EXTRACT_UTA_EXON_GAP_INFO.out.uta_exon_gap_info_csv)
    emit:
        gff_and_uta_exon_gap_info = COMBINE_GFF_AND_UTA_EXON_GAP_INFO.out.gff_and_uta_exon_gap_info
}