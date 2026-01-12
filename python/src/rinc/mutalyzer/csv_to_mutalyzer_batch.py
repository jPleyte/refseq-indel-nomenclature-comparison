'''
Read a list of variants and write them out in format suitable for submission to Mutilizer''s Batch Processor

Created on Jan 11, 2026

@author: pleyte
'''
import argparse
import logging.config
from rinc.util.log_config import LogConfig
from rinc.variant_transcript import VariantTranscript
import csv
from pathlib import Path

class MutalyzerBatch(object):
    '''
    Utility class for extracting variants from the the gap-variant csv and writing out a list of g. variants that can be submitted to Mutilizer's Batch Processor  
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        
    def get_variants(self, variants_file: str):
        """
        Read the csv file that has the variants that will be processed. 
        """
        variants = [] 
        with open(variants_file, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                variants.append(VariantTranscript(row['chromosome'], 
                                                  int(row['position']),
                                                  row['reference'], 
                                                  row['alt'],
                                                  row['cdna_transcript'], 
                                                  g_dot=row['g_dot']))
        
        self._logger.debug(f"Read {len(variants)} variants from {variants_file}")
        return variants
    
    def write_transcripts(self, out_directory, variants: list[VariantTranscript], batch_size=50):
        """
        Write variant transcripts out to batch files in the out_directory. Each file will have max batch_size lines
        """
        batch_n = -1
        
        batches = [variants[i:i + batch_size] for i in range(0, len(variants), batch_size)]
        
        for batch in batches:
            batch_n= batch_n + 1
            directory = Path(out_directory)
            file = f"mutalyzer_batch_{batch_n}.txt"
            full_path = directory / file
            
            with open(full_path, "w") as f:
                for v in batch:
                    # f.write(f"{v.g_dot}\t{v.cdna_transcript}\n")
                    refseq_chromosome, g_dot_only = v.g_dot.split(':')
                    f.write(f"{refseq_chromosome}:{g_dot_only}\n")
                
                self._logger.info(f"Wrote {len(batch)} variant transcripts to {full_path}")
        

def _parse_args():
    parser = argparse.ArgumentParser(description='Read a list of variants and write them out in format suitable for submission to Mutilizer''s Batch Processor')

    parser.add_argument('--variants', help='Variants (csv)', required=True)
    parser.add_argument('--out_dir', help='Directory to write batch files', required=True)

    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    mb = MutalyzerBatch()
    variants = mb.get_variants(args.variants)
    
    mb.write_transcripts(args.out_dir, variants)
    
if __name__ == '__main__':
    main()