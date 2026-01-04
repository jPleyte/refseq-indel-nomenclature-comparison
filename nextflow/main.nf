#!/usr/bin/env nextflow

// Import the process from your modules folder
include { findGaps } from './modules/local/find_gaps.nf'

workflow {
    main:
    findGaps()
}