#!/usr/bin/env nextflow

/*
 * Instead of installing Mutalyzer locally or spending loads of time performing api requests, mutliizer annoations
 * are stored in a local csv file before the workflow is launched. This process looks up each transcript in the local 
 * db and writes out a csv with the ones it finds.   
 */

 process writeMutalyzerNomenclatureToCsv {
     publishDir "${params.outdir}", mode: 'symlink'
 
     input:
     path variants_csv
     path mutalyzer_cache
 
     output:
     path "mutalyzer_nomenclature.csv", emit: mutalyzer_nomenclature
 
     script:
     """
     python -m rinc.mutalyzer.mutalyzer_nomenclature \
         --variants ${variants_csv} \
         --mutalyzer_cache ${mutalyzer_cache} \
         --out mutalyzer_nomenclature.csv
     """
 }