#!/usr/bin/env nextflow

/*
 * Run SnpEff on the provided snpeff input file.
 */
process runSnpEff {
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path variants_vcf
    
    output:
    path "snpeff_output.vcf", emit: snpeff_vcf
    path "snpeff_nomenclature.tsv", emit: snpeff_tsv

    script:
    """ 
    # snpeff options: -verbose -noStats -csvStats x.csv -canon
    java -Xmx8g -jar ${params.snpeff_jar} \
         ann \
         -noStats \
         GRCh37.p13 \
         ${variants_vcf} > snpeff_output.vcf

    cat snpeff_output.vcf | perl ${params.snpeff_oneperline_pl} | java -jar ${params.snpsift_jar} extractFields \
         - \
         -e "." \
         CHROM POS ID REF ALT "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].IMPACT" "ANN[*].BIOTYPE" > snpeff_nomenclature.tsv
    """
}
