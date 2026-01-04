#!/usr/bin/env nextflow

// Import the process from your modules folder
include { findGapVariants } from './modules/local/find_gap_variants.nf'

workflow {
    main:
    uta_schema = channel.value(params.uta_schema)
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
    
    findGapVariants(uta_schema, fasta_ch)
}