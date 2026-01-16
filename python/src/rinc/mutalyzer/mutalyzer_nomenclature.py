'''
Annotate variant_helper with results from csvs downloaded from Mutalyzer Batch Processor

Created on Jan 11, 2026

@author: pleyte
'''

import logging.config
import argparse
from rinc.util.log_config import LogConfig
from rinc.mutalyzer.mutalyzer_cache import MytalyzerCache
from rinc.variant_transcript import VariantTranscript
import csv
from rinc.util.pdot import PDot
from rinc.io import variant_helper
import pandas as pd

class MutalyzerNomenclature(object):
    '''
    Annotate variant_helper with results from Mutalyzer Batch Processor
    '''
    def __init__(self, mutalyzer_cache_file):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._mutalyzer_cache = MytalyzerCache(mutalyzer_cache_file)
        self._pDot = PDot()
    
    def _get_full_g_dot(self, mutalyzer_result):
        """
        """
        refseq_transcript = str(mutalyzer_result.refseq_chromosome) 
        g_dot = str(mutalyzer_result.g_dot)
        return f"{refseq_transcript}:{g_dot}"
        
    def get_annotated_variants(self, variant_transcripts: list[VariantTranscript]):
        annotated_variants = []
        
        for v in variant_transcripts:
            mutalyzer_result = self._get_result(v)
            if not mutalyzer_result.empty:
                v.c_dot = None if pd.isna(mutalyzer_result.c_dot) else str(mutalyzer_result.c_dot) 
                v.g_dot = None if pd.isna(mutalyzer_result.g_dot) else self._get_full_g_dot(mutalyzer_result)
                v.p_dot3 = None if pd.isna(mutalyzer_result.p_dot) else str(mutalyzer_result.p_dot)                
                v.protein_transcript = None if pd.isna(mutalyzer_result.protein_transcript) else str(mutalyzer_result.protein_transcript)
                
                if v.p_dot3:
                    v.p_dot1 = self._pDot.get_p_dot1(v.protein_transcript, v.p_dot3)

                v.additional_fields['mut.status'] = str(mutalyzer_result.status)
                annotated_variants.append(v)
        
        self._logger.info(f"Found {len(annotated_variants)} of {len(variant_transcripts)} in Mutalyzer cache")
        return annotated_variants
            
    def _get_result(self, variant: VariantTranscript) -> pd.DataFrame:
        """
        Look for the variant transcript in the mutalyzer cache
        """
        if not variant.g_dot:
            raise ValueError("Mutalyzer requires a g_dot in the variant source")
        
        refseq_chromosome, g_dot = variant.g_dot.split(':')

        return self._mutalyzer_cache.find(refseq_chromosome, g_dot, variant.cdna_transcript)

    def write(self, output_file: str, variant_transcripts: list[VariantTranscript]):
        """
        """
        headers = ['chromosome', 'position', 'reference', 'alt',
                   'cdna_transcript', 'mut.protein_transcript', 
                   'mut.g_dot', 'mut.c_dot', 'mut.p_dot1', 'mut.p_dot3',
                   'mut.status']
        
        with open(output_file, 'w', newline='') as output:
            writer = csv.writer(output)
            writer.writerow(headers)

            for v in variant_transcripts:
                writer.writerow([v.chromosome, v.position, v.reference, v.alt, 
                                 v.cdna_transcript, v.protein_transcript,
                                 v.g_dot, v.c_dot, v.p_dot1, v.p_dot3,
                                 v.additional_fields['mut.status']])
            
        self._logger.info(f"Wrote {len(variant_transcripts)} variant transcripts to {output_file}")        

def _parse_args():
    parser = argparse.ArgumentParser(description='Lookup variant_helper in local mutalyzer db and write annotations to output file')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--variants", help="File with variant_helper (csv)", required=True)
    parser.add_argument("--mutalyzer_cache", help="File with variant (csv)", required=True)
    parser.add_argument("--out", help="output file (csv)", dest="output", required=True)
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()
    
    mn = MutalyzerNomenclature(args.mutalyzer_cache)
    
    # Read variant list from file
    variants = variant_helper.get_variants(args.variants)
    logging.debug(f"Read {len(variants)} variants from {args.variants}")
    
    # Lookup variant in mutalyzer db, update nomenclature, and return only the ones that are found
    variants_with_annoation = mn.get_annotated_variants(variants)
    
    mn.write(args.output, variants_with_annoation)
    

if __name__ == '__main__':
    main()