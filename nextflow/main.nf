#!/usr/bin/env nextflow

// Import the process from your modules folder
include { findGaps } from './modules/local/find_gaps.nf'

workflow {
    main:
    // params.fasta should be defined in nextflow.config or on the command line
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)

    findGaps(fasta_ch)
}