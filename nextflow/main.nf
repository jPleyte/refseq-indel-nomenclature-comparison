#!/usr/bin/env nextflow

// Import the processes from your modules folder
include { validateParameters; paramsSummaryLog } from 'plugin/nf-validation'
include { findGapVariants } from './modules/local/variants/find_gap_variants.nf'
include { getTfxVariants } from './modules/local/variants/get_tfx_variants.nf'
include { writeHgvsNomenclatureToCsv} from './modules/local/write_hgvs_nomenclature_to_csv.nf'
include { csvToAvinput } from './modules/local/csv_to_avinput.nf'
include { runAnnovar } from './modules/local/run_annovar.nf'
include { writeAnnovarNomenclatureToCsv } from './modules/local/write_annovar_nomenclature_to_csv.nf'
include { csvToVcf } from './modules/local/csv_to_vcf.nf'
include { runSnpEff } from './modules/local/run_snpeff.nf'
include { writeSnpEffNomenclatureToCsv } from './modules/local/write_snpeff_nomenclature_to_csv.nf'
include { runVep as runVepRefseq } from './modules/local/run_vep.nf'
include { runVep as runVepHg19 } from './modules/local/run_vep.nf'
include { writeVepNomenclatureToCsv as writeVepRefseqNomenclatureToCsv } from './modules/local/write_vep_nomenclature_to_csv.nf'
include { writeVepNomenclatureToCsv as writeVepHg19NomenclatureToCsv } from './modules/local/write_vep_nomenclature_to_csv.nf'
include { writeTfxNomenclatureToCsv } from './modules/local/write_tfx_nomenclature_to_csv.nf'
include { writeCgdNomenclatureToCsv } from './modules/local/write_cgd_nomenclature_to_csv.nf'
include { joinAndCompare } from './modules/local/join_and_compare.nf'
include { performAnalysis } from './modules/local/perform_analysis.nf'

workflow {
    main:
    uta_schema = channel.value(params.uta_schema)
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)

    def vep_sequence_modes = [
        refseq: "--use_transcript_ref",
        reference: "--use_given_ref"
    ]

    validateParameters()

    println "Using variant source: ${params.variant_source}"
    if (params.variant_source != 'gap_query') {
        println "Using variant source file: ${params.variant_source_file}"
    } 

    // Parameter validation
    if (params.variant_source_file != null && params.variant_source == 'gap_query') {
        error "You provided an input file but chose gap_query as the input method, which does not use a file."
    } else if (params.variant_source_file == null && (params.variant_source == 'csv' || params.variant_source == 'tfx')) {
        error "You indicated ${params.variant_source} as the variant_source, but did not provide a variant_source_file"
    }

    // Read variants 
    def ch_variants
    if (params.variant_source == 'csv') {
        ch_variants = params.variant_source_file
    } 
    else if (params.variant_source == 'gap_query') {
        ch_variants = findGapVariants(uta_schema, fasta_ch)
    }
    else if (params.variant_source == 'tfx') {
        ch_variants = getTfxVariants(params.variant_source_file)
    }
    else {
        error "Unknown variant source: ${params.variant_source}"
    }

    // Generate nomenclature using hgvs package. 
    // Do not generate hgvs nomenclature when input source is 'tfx' but tfx uses the same code.
    def hgvs_nomenclature = channel.empty()
     if (params.variant_source != 'tfx') {
        hgvs_nomenclature = writeHgvsNomenclatureToCsv(fasta_ch, ch_variants)
     }
    

    // Convert variant list to annovar avinput file 
    csvToAvinput(ch_variants)

    // Run annovar on the avinput file
    runAnnovar(csvToAvinput.out.annovar_avinput)

    // Extract annovar nomenclature and write to new csv
    writeAnnovarNomenclatureToCsv(runAnnovar.out.multianno)
    
    // Convert variant list to vcf to be used by SnpEff and VEP
    csvToVcf(ch_variants)

    // Run SnpEfff on vcf 
    runSnpEff(csvToVcf.out.vcf)

    // Extract SnpEff nomenclature and write to new csv
    writeSnpEffNomenclatureToCsv(runSnpEff.out.snpeff_tsv)

    // Run VEP onece using coding sequence for reference andusing hg19 for reference 
    runVepRefseq(csvToVcf.out.vcf, params.vep_fasta, 'refseq', vep_sequence_modes.refseq)
    runVepHg19(csvToVcf.out.vcf, params.vep_fasta, 'hg19', vep_sequence_modes.reference)

    // extract VEP nomenclature and write to new csv
    writeVepRefseqNomenclatureToCsv(runVepRefseq.out.vep_output, 'refseq')
    writeVepHg19NomenclatureToCsv(runVepHg19.out.vep_output, 'hg19')

    def tfx_nomenclature = channel.empty()
    if (params.variant_source == 'tfx') {
        tfx_nomenclature = writeTfxNomenclatureToCsv(fasta_ch, params.variant_source_file)
    }

    def cgd_nomenclature = channel.empty()
    if (params.cgd_export_csv != null) {
        cgd_nomenclature = writeCgdNomenclatureToCsv(ch_variants)
    }

    // Compare hgvs and annovar, join hgvs, annovar, and gaps file into final output
    joinAndCompare(ch_variants, 
                   hgvs_nomenclature.ifEmpty([]),
                   writeAnnovarNomenclatureToCsv.out.annovar_nomenclature,
                   writeSnpEffNomenclatureToCsv.out.snpeff_nomenclature,
                   writeVepRefseqNomenclatureToCsv.out.vep_nomenclature,
                   writeVepHg19NomenclatureToCsv.out.vep_nomenclature,
                   tfx_nomenclature.ifEmpty([]),
                   cgd_nomenclature.ifEmpty([]))

    performAnalysis(hgvs_nomenclature.ifEmpty([]),
                    writeAnnovarNomenclatureToCsv.out.annovar_nomenclature,
                    writeSnpEffNomenclatureToCsv.out.snpeff_nomenclature,
                    writeVepRefseqNomenclatureToCsv.out.vep_nomenclature,
                    writeVepHg19NomenclatureToCsv.out.vep_nomenclature,
                    tfx_nomenclature.ifEmpty([]),
                    cgd_nomenclature.ifEmpty([]))
}