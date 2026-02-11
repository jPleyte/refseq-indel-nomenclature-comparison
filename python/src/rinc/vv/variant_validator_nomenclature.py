'''
After submitting a list of variants to VariantValidator.org we get back a tsv.  
This module takes the received tsv file extracts the relevant nomenclature
and writes it out in our common csv format. 

Created on Jan 30, 2026

@author: pleyte
'''
import argparse
import logging.config
from rinc.util.log_config import LogConfig
import pandas as pd
from rinc.util.pdot import PDot
from collections import Counter

class VariantValidatorNomenclature(object):
    '''
    Convert VariantValidator data to common csv format so it can be included in pipeline analysis.  
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._p_dot_converter = PDot()
        self._nomenclature_counter = Counter()
    
    def _get_cdna_transcript_and_c_dot(self, hgvs_transcript: str):
        """
        Extract transcript and c. from the variant validator HGVS_transcript field         
        """
        if type(hgvs_transcript) == float:
            print("jdebug")
        elif not hgvs_transcript:
            # Field will be blank when our transcript filter doesn't match anything
            return None, None
             
        cdna_transcript, c_dot = hgvs_transcript.split(':')
        if not cdna_transcript or not c_dot:
            raise ValueError("one or more cdna fields are blank")
        return cdna_transcript, c_dot
        
    def _get_protein_transcript_and_p_dot(self, hgvs_predicted_protein):
        """
        Extract protein transcript and p. from the variant validator HGVS_Predicted_Protein field
        """
        if ':' not in hgvs_predicted_protein:
            print("jDebug")
            
        protein_transcript, p_dot = hgvs_predicted_protein.split(':')
        if not protein_transcript or not p_dot:
            raise ValueError("one or more protein fields are blank")
        return protein_transcript, p_dot
        
    def _get_variant_transcript(self, row):
        """
        Extract values from VariantValidator row 
        """
        if 'None of the specified transcripts' in row.Warnings:
            self._nomenclature_counter['nothing_due_to_transcript_filter'] += 1            
            cdna_transcript = c_dot = protein_transcript = p_dot1 = p_dot3 = None
        else:
            cdna_transcript, c_dot = self._get_cdna_transcript_and_c_dot(row.HGVS_transcript)
            protein_transcript, p_dot3_raw = self._get_protein_transcript_and_p_dot(row.HGVS_Predicted_Protein)
            p_dot3 = self._p_dot_converter.get_remove_parenthesis(protein_transcript, p_dot3_raw, True)        
            p_dot1 = self._p_dot_converter.get_p_dot1(protein_transcript, p_dot3)
            self._nomenclature_counter['transcript_nomenclature_found'] += 1
        
        return {
            'chromosome': row.GRCh37_CHR,
            'position':  row.GRCh37_POS,
            'reference':  row.GRCh37_REF,
            'alt':  row.GRCh37_ALT,
            'gene':  row.Gene_Symbol,
            'g_dot':  row.HGVS_Genomic_GRCh37,
            'cdna_transcript': cdna_transcript, 
            'c_dot': c_dot,
            'protein_transcript': protein_transcript,
            'p_dot1': p_dot1,
            'p_dot3': p_dot3            
        }
        
    def get_variant_nomenclature(self, vv_file: str):
        df = pd.read_csv(vv_file, sep='\t', comment='#', dtype=str, keep_default_na=False)
        self._logger.info(f"Read {df.shape[0]} rows from {vv_file}")
        variant_transcripts = []
        for row in df.itertuples():
            variant_transcript = self._get_variant_transcript(row)
            variant_transcripts.append(variant_transcript)
        
        self._logger.info(f"Nomenclature counter: {self._nomenclature_counter}")
        return variant_transcripts
    
    def write(self, output_file, variant_transcripts: list[dict]):
        output_df = pd.DataFrame(variant_transcripts)
        output_df.to_csv(output_file, index=False)
        self._logger.info(f"Wrote {output_df.shape[0]} rows to {output_file}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Read VariantValidator file and write out new csv')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--vv_input", help="Variant Validator batch received by email (tsv)", required=True)
    parser.add_argument("--nomenclature_output", help="Nomenclature output file (csv)", required=True)
    
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    vvn = VariantValidatorNomenclature()
    
    variant_transcripts = vvn.get_variant_nomenclature(args.vv_input)
    vvn.write(args.nomenclature_output, variant_transcripts)

if __name__ == '__main__':
    main()