#!/usr/bin/env nextflow

/*
 * Filter a csv with variants and nomenclature and keep only the variant-transcripts that match what is in the filter_file
 */
process filterNomenclature {
    tag "${input_csv.baseName}"
    label 'process_low'

    input:
    path input_csv
    path filter_file // Staged as optional via the workflow logic below

    output:
    path "${input_csv.baseName}_filtered.csv", emit: filtered_csv

    script:
    if (filter_file.name != 'NO_FILTER') 
        """
        python -m rinc.io.filter_nomenclature \
            --input ${input_csv} \
            --filter ${filter_file} \
            --output ${input_csv.baseName}_filtered.csv
        """
    else
        """
        cp ${input_csv} ${input_csv.baseName}_filtered.csv
        """
}