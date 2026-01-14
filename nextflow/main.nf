#!/usr/bin/env nextflow

// Import the processes from your modules folder
include { findGapVariants } from './modules/local/find_gap_variants.nf'
include { writeHgvsNomenclatureToCsv} from './modules/local/write_hgvs_nomenclature_to_csv.nf'
include { csvToAvinput } from './modules/local/csv_to_avinput.nf'
include { runAnnovar } from './modules/local/run_annovar.nf'
include { writeAnnovarNomenclatureToCsv } from './modules/local/write_annovar_nomenclature_to_csv.nf'
include { csvToVcf } from './modules/local/csv_to_vcf.nf'
include { runSnpEff } from './modules/local/run_snpeff.nf'
include { writeSnpEffNomenclatureToCsv } from './modules/local/write_snpeff_nomenclature_to_csv.nf'
include { writeMutalyzerNomenclatureToCsv } from './modules/local/write_mutalyzer_nomenclature_to_csv.nf'
include { runVep as runVepRefseq } from './modules/local/run_vep.nf'
include { runVep as runVepHg19 } from './modules/local/run_vep.nf'
include { writeVepNomenclatureToCsv as writeVepRefseqNomenclatureToCsv } from './modules/local/write_vep_nomenclature_to_csv.nf'
include { writeVepNomenclatureToCsv as writeVepHg19NomenclatureToCsv } from './modules/local/write_vep_nomenclature_to_csv.nf'
include {find_missing_mutalyzer_variants } from './modules/local/find_missing_mutalyzer_variants.nf'


include { joinAndCompare } from './modules/local/join_and_compare.nf'

workflow {
    main:
    uta_schema = channel.value(params.uta_schema)
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
    
    def vep_sequence_modes = [
        refseq: "--use_transcript_ref",
        reference: "--use_given_ref"
    ]

    // Create a variant list 
    findGapVariants(uta_schema, fasta_ch)

    // Generate nomenclature using hgvs package
    writeHgvsNomenclatureToCsv(findGapVariants.out.gaps_and_variants)

    // Convert variant list to annovar avinput file 
    csvToAvinput(findGapVariants.out.gaps_and_variants)

    // Run annovar on the avinput file
    runAnnovar(csvToAvinput.out.annovar_avinput)

    // Extract annovar nomenclature and write to new csv
    writeAnnovarNomenclatureToCsv(runAnnovar.out.multianno)
    
    // Convert variant list to vcf to be used by SnpEff and VEP
    csvToVcf(findGapVariants.out.gaps_and_variants)

    // Run SnpEfff on vcf 
    runSnpEff(csvToVcf.out.vcf)

    // Extract SnpEff nomenclature and write to new csv
    writeSnpEffNomenclatureToCsv(runSnpEff.out.snpeff_tsv)

    // Lookup Muttalizaer annotatinos in local cache
    writeMutalyzerNomenclatureToCsv(findGapVariants.out.gaps_and_variants, params.mutilizer_cache)

    // Run VEP onece using coding sequence for reference andusing hg19 for reference 
    runVepRefseq(csvToVcf.out.vcf, params.vep_fasta, 'refseq', vep_sequence_modes.refseq)
    runVepHg19(csvToVcf.out.vcf, params.vep_fasta, 'hg19', vep_sequence_modes.reference)

    // extract VEP nomenclature and write to new csv
    writeVepRefseqNomenclatureToCsv(runVepRefseq.out.vep_output, 'refseq')
    writeVepHg19NomenclatureToCsv(runVepHg19.out.vep_output, 'hg19')

    // Compare hgvs and annovar, join hgvs, annovar, and gaps file into final output
    joinAndCompare(findGapVariants.out.gaps_and_variants, 
                   writeHgvsNomenclatureToCsv.out.hgvs_nomenclature,
                   writeAnnovarNomenclatureToCsv.out.annovar_nomenclature,
                   writeSnpEffNomenclatureToCsv.out.snpeff_nomenclature,
                   writeMutalyzerNomenclatureToCsv.out.mutalyzer_nomenclature,
                   writeVepRefseqNomenclatureToCsv.out.vep_nomenclature,
                   writeVepHg19NomenclatureToCsv.out.vep_nomenclature)

    find_missing_mutalyzer_variants(joinAndCompare.out.gap_nomenclature_all_transcripts_all_fields, params.mutilizer_cache)
}