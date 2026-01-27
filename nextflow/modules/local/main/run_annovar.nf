#!/usr/bin/env nextflow

/*
 * Run ANNOVAR on the provided avinput file.
 */
process runAnnovar {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path annovar_avinput
    
    output:
    path "annovar.hg19_multianno.txt", emit: multianno
    
    script:

    """    
    \$ANNOVAR_HOME/table_annovar.pl ${annovar_avinput} \
    \$ANNOVAR_HOME/humandb/ \
    --buildver hg19 \
    --out annovar \
    --protocol refGeneWithVer,ccdsGene \
    --operation g,g \
    --nastring . \
    --polish \
    --remove \
    --argument '--splicing_threshold 5 --exonicsplicing --transcript_function --separate,--splicing_threshold 5 --exonicsplicing --transcript_function --separate'
    """
}
