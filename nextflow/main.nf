#!/usr/bin/env nextflow

// Import the processes from your modules folder
include { findGapVariants } from './modules/local/find_gap_variants.nf'
include { hgvsNomenclature} from './modules/local/hgvs_nomenclature.nf'
include { csvToAvinput } from './modules/local/csv_to_avinput.nf'
include { runAnnovar } from './modules/local/run_annovar.nf'
include { convertAnnovarMultiannoToCsv } from './modules/local/convert_annovar_multianno_to_csv.nf'
include { csvToVcf } from './modules/local/csv_to_vcf.nf'
include { runSnpEff } from './modules/local/run_snpeff.nf'
include { convertSnpEffToCsv } from './modules/local/convert_snpeff_to_csv.nf'

//include { csvToVep } from './modules/local/csv_to_vep.nf'
//include { runVep } from './modules/local/run_vep.nf'


include { joinAndCompare } from './modules/local/join_and_compare.nf'

workflow {
    main:
    uta_schema = channel.value(params.uta_schema)
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
    
    // Create a variant list 
    findGapVariants(uta_schema, fasta_ch)

    // Generate nomenclature using hgvs package
    hgvsNomenclature(findGapVariants.out.gaps_and_variants)

    // Convert variant list to annovar avinput file 
    csvToAvinput(findGapVariants.out.gaps_and_variants)

    // Run annovar on the avinput file
    runAnnovar(csvToAvinput.out.annovar_avinput)

    // Extract annovar nomenclature and write to new csv
    convertAnnovarMultiannoToCsv(runAnnovar.out.multianno)
    
    // Convert variant list to vcf 
    csvToVcf(findGapVariants.out.gaps_and_variants)

    // Run SnpEfff on vcf 
    runSnpEff(csvToVcf.out.vcf)

    // Extract SnpEff nomenclature and write to new csv
    convertSnpEffToCsv(runSnpEff.out.snpeff_tsv)

    // Convert variant list to VEP input file
    // csvToVep(findGapVariants.out.gaps_and_variants)

    // Compare hgvs and annovar, join hgvs, annovar, and gaps file into final output
    joinAndCompare(findGapVariants.out.gaps_and_variants, 
                   hgvsNomenclature.out.hgvs_nomenclature,
                   convertAnnovarMultiannoToCsv.out.annovar_nomenclature,
                   convertSnpEffToCsv.out.snpeff_nomenclature)
}