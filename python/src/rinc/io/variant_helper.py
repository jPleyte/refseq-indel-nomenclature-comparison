'''
Created on Jan 14, 2026

@author: pleyte
'''
from rinc.variant_transcript import VariantTranscript
import csv
import logging


def get_variants(variants_file: str):
    """
    Read the csv file that has the variants that will be processed. 
    """
    variants = [] 
    with open(variants_file, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            v = VariantTranscript(row['chromosome'],
                                  int(row['position']),
                                  row['reference'],
                                  row['alt'],
                                  row['cdna_transcript'])
            if 'g_dot' in row:
                v.g_dot = row['g_dot']
                              
            variants.append(v)

    return variants


def write_variants(out_filename: str, variants: list[VariantTranscript]):
    """
    Write variants to csv file
    """
    headers = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript', 'g_dot']
    
    with open(out_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        
        for v in variants: 
            writer.writerow([v.chromosome, v.position, v.reference, v.alt, v.cdna_transcript, v.g_dot])

        
def write_variant_transcripts(out_filename: str, variants: list[VariantTranscript], additional_fields=[], field_suffix=""):
    """
    Write a list of VariantTranscript to csv 
    Each field label will have the field_suffix appended to it
    """
    
    key_headers = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript' ] 
    nomenclature_headers = ['c_dot', 'exon', 'g_dot', 'gene', 'genomic_region_type', 'p_dot1', 'p_dot3', 'protein_transcript', 'protein_variant_type']
    nomenclature_headers.extend(additional_fields)
    suffixed_headers = [x + "." + field_suffix for x in nomenclature_headers]
    all_headers = key_headers + suffixed_headers 
    
    rows = 0
    with open(out_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(all_headers)
        
        for v in variants:
            rows += 1 
            row = [v.chromosome,
                   v.position,
                   v.reference,
                   v.alt,
                   v.cdna_transcript,
                   v.c_dot,
                   v.exon,
                   v.g_dot,
                   v.gene,
                   v.genomic_region_type,
                   v.p_dot1,
                   v.p_dot3,
                   v.protein_transcript,
                   v.protein_variant_type]
            
            for x in additional_fields:
                row.append(v.additional_fields[x])
                
            writer.writerow(row)
        
    logging.info(f"Wrote {rows} to {out_filename}")
