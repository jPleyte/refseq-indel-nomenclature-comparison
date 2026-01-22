'''
Read a CGD export (csv) and write out a distinct list of variant and transcripts
Created on Jan 21, 2026

@author: pleyte
'''
import argparse
import logging.config
import pandas as pd
from rinc.util.log_config import LogConfig


class CgdCsvToVariantTranscriptCsv(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)

    def get_variants(self, input_filename):
        """
        Read a csv containing variants and transcript nomenclature 
        """
        df = pd.read_csv(input_filename, usecols=['chromosome', 'position_start', 'reference_base', 'variant_base', 'cdna_transcript']).drop_duplicates()
        
        # Remove "chr" prefix from chromosomes
        df['chromosome'] = df['chromosome'].astype(str).str.replace('^chr', '', case=False, regex=True)
        
        self._logger.info(f"Read {df.shape[0]} rows from {input_filename}")
        return df
    
    def write(self, output_filename, variants_df: pd.DataFrame):
        """
        """
        variants_df.rename(columns={'chromosome': 'chromosome',
                                    'position_start': 'position',
                                    'reference_base': 'reference',
                                    'variant_base': 'alt',
                                    'cdna_transcript': 'cdna_transcript'}, inplace=True)
        
        
        variants_df.to_csv(output_filename, index=False)
        self._logger.info(f"Wrote {variants_df.shape[0]} rows to {output_filename}")


def _parse_args():
    parser = argparse.ArgumentParser(description='Convert sql dump to variant transcript list')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--input", help="SQL dump (csv)", required=True)    
    parser.add_argument("--output", help="Variant transcript list (csv)", required=True)
    
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()
    
    cgd_to_csv = CgdCsvToVariantTranscriptCsv()
    variant_transcripts = cgd_to_csv.get_variants(args.input)
    cgd_to_csv.write(args.output, variant_transcripts)
    
if __name__ == '__main__':
    main()        