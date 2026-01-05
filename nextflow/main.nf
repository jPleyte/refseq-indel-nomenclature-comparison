#!/usr/bin/env nextflow

// Import the processes from your modules folder
include { findGapVariants } from './modules/local/find_gap_variants.nf'
include { hgvsNomenclature} from './modules/local/hgvs_nomenclature.nf'
include { csvToAvinput } from './modules/local/csv_to_avinput.nf'
include { runAnnovar } from './modules/local/run_annovar.nf'

workflow {
    main:
    uta_schema = channel.value(params.uta_schema)
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
    
    // Create a variant list 
    findGapVariants(uta_schema, fasta_ch)

    // Generate nomenclature using hgvs package
    // hgvsNomenclature(findGapVariants.out.gaps_and_variants)

    // Convert variant list to annovar avinput file 
    csvToAvinput(findGapVariants.out.gaps_and_variants)

    runAnnovar(csvToAvinput.out.annovar_avinput)
}