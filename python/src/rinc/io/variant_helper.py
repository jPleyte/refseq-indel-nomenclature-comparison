'''
Created on Jan 14, 2026

@author: pleyte
'''
from rinc.variant_transcript import VariantTranscript
import csv

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