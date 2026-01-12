#!/usr/bin/env nextflow

/*
 * Instead of installing Mutilizer locally or spending loads of time performing api requests, mutliizer annoations
 * are stored in a local csv file before the workflow is launched. This process looks up each transcript in the local 
 * db and writes out a csv with the ones it finds.   
 */

 process addMutalizerAnnoation {
     publishDir "${params.outdir}", mode: 'symlink'
 
     input:
     path variants_csv
     path mutilizer_cache
 
     output:
     path "mutalizer_nomenclature.csv", emit: mutalizer_nomenclature
 
     script:
     """
     python -m rinc.mutalyzer.mutalyzer_nomenclature \
         --variants ${variants_csv} \
         --mutalyzer_cache ${mutilizer_cache} \
         --out mutalizer_nomenclature.csv
     """
 }