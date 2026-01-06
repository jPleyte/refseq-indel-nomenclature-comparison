#!/usr/bin/env nextflow

/*
 * Run ANNOVAR on the provided avinput file.
 */
process runAnnovar {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path annovar_avinput
    
    output:
    path "annovar.hg19_multianno.csv", emit: multianno
    
    script:
    """    
    \$ANNOVAR_HOME/table_annovar.pl ${annovar_avinput} \
    \$ANNOVAR_HOME/humandb/ \
    --buildver hg19 \
    --out annovar \
    --protocol refGeneWithVer \
    --operation g \
    --nastring . \
    --polish \
    --csvout \
    --remove \
    --argument '--splicing_threshold 5 --exonicsplicing --transcript_function --separate'
    """
}
