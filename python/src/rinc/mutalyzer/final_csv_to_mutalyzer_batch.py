'''
Read the final output of the workflow and pick out variants that we want to submit to to the 
Mutalyzer Batch Processor. 

I don't want to install Mutalyzer locally but I am going to have a lot of variants to look at which
would make sending individual REST requests to Mutalyzer take a long time. So what i'm doing is 
scanning the output of the workflow for variants that meet certain criteria and then putting 
them in batch files that i can submit to the Mutalyzer Batch Processor. I'll then take the 
results and impomrt them into my local Mutalyzer cache (using mutalyzer_cache.py) and 
re-run the workflow which will then add Mutalyzer nomenclature for those variants to the final output.   

Created on Jan 11, 2026

@author: pleyte
'''
import argparse
import logging.config
from rinc.util.log_config import LogConfig
from rinc.variant_transcript import VariantTranscript
import csv
from pathlib import Path
from rinc.mutalyzer.mutalyzer_cache import MytalyzerCache
from collections import Counter
from rinc.util import chromosome_map

class FinalCsvToMutalyzerBatch(object):
    '''
    Utility class for extracting variants from the the gap-variant csv and writing out a list of g. variants that can be submitted to Mutalyzer's Batch Processor  
    '''

    def __init__(self, mutalyzer_cache_file=None):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._mutalyzer_cache = self._load_mutalyzer_cache(mutalyzer_cache_file)
        self._criteria_counter = Counter()
        
    def _load_mutalyzer_cache(self, mutalyzer_cache_file: str):
        """
        Load the Mutalyzer cache. This is optional. It allows us to check the local cache for results before 
        adding them to the batch so we don't request annotation for variants we've already requested. 
        """
        if not mutalyzer_cache_file:
            return None
        
        return MytalyzerCache(mutalyzer_cache_file)
    
    def _is_result_already_known(self, row: dict) -> bool:
        """
        Check the Mutalyzer cache for the variant and return true if it is found. 
        Returns False for all variants if the Mutalyzer cache does not exist. 
        """
        if not self._mutalyzer_cache:
            return False
        
        refseq_chromosome_a = chromosome_map.get_refseq(row['chromosome'])
        refseq_chromosome_b = row['alt_ac']
        refseq_chromosome_c, g_dot = row['g_dot'].split(':')
        assert refseq_chromosome_a == refseq_chromosome_b == refseq_chromosome_c, f"Chromosomes do not match for {row['chromosome']}-{row['position']}-{row['reference']}-{row['alt']}-{row['cdna_transcript']}: {refseq_chromosome_a} != {refseq_chromosome_b} != {refseq_chromosome_c}"
        
        refseq_cdna_transcript = row['cdna_transcript']
        assert refseq_cdna_transcript, f"cnda transcript is blank for {row['chromosome']}-{row['position']}-{row['reference']}-{row['alt']}"
            
        # Check the db for this variant and transcript. The first serach uses the transcript. 
        # But If mutalyzer doesn't have any info for this variant then it will  return an empty result which is stored in the local db without a transcript so
        # a second search is performed without the transcript.   
        if self._mutalyzer_cache.find(refseq_chromosome_a, g_dot, refseq_cdna_transcript):
            return True
        elif self._mutalyzer_cache.find(refseq_chromosome_a, g_dot, None):
            
            return True
        else:
            return False
        
    def _meets_criteria(self, row: dict) -> bool:
        """
        Look at the row and make sure it has c. from hgvs/uta and annovar. As long as
        it has that much then it's worth submitting to Mutalyzer and including in 
        our results. 
        """
        
        
        if row['hu.c_dot'] and row['annovar.c_dot'] and self._is_result_already_known(row):
            self._logger.debug(f"Skipping already known {row['chromosome']}-{row['position']}-{row['reference']}-{row['alt']}-{row['cdna_transcript']}")
            self._criteria_counter['already_cached'] += 1
            return False
        elif row['hu.c_dot'] and row['annovar.c_dot']:
            self._logger.debug(f"Keeping {row['chromosome']}-{row['position']}-{row['reference']}-{row['alt']}-{row['cdna_transcript']} with {row['hu.c_dot']} and {row['annovar.c_dot']}")
            self._criteria_counter['kept'] += 1
            return True         
        else:
            self._logger.debug(f"Skipping {row['chromosome']}-{row['position']}-{row['reference']}-{row['alt']}-{row['cdna_transcript']} with {row['hu.c_dot']} and {row['annovar.c_dot']}")
            self._criteria_counter['lacking_sufficient_c_dot'] += 1
            return False
        
    def get_variants(self, variants_file: str):
        """
        Read the csv file that has the variants that will be processed. 
        """
        result_counter = Counter()
        variants_to_request = set()
        
        with open(variants_file, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                result_counter['total'] += 1
                if row['position'] == '46694318' and row['reference'] == 'A' and row['alt'] == 'G':
                    print("jDebug")
                    
                if self._meets_criteria(row):                    
                    variant_string = self._get_mutalyzer_batch_variant(VariantTranscript(row['chromosome'], 
                                                                                         int(row['position']),
                                                                                         row['reference'], 
                                                                                         row['alt'],
                                                                                         row['cdna_transcript'], 
                                                                                         g_dot=row['g_dot']))
                    if variant_string in variants_to_request:
                        result_counter['duplicate'] += 1
                    else:
                        result_counter['unique'] += 1
                        variants_to_request.add(variant_string)
                
                
        self._logger.info(f"Read {result_counter['total']} variants from {variants_file}")
        self._logger.info(f"Criteria: {self._criteria_counter}")
        self._logger.info(f"Kept {result_counter['unique']}, discarded {result_counter['duplicate']} duplicates")
        
        
        return list(variants_to_request)
    
    def _get_mutalyzer_batch_variant(self, v: VariantTranscript) -> str:
        """
        Conver the VariantTranscript to the g. string that will be put in the batch file
        """
        refseq_chromosome, g_dot_only = v.g_dot.split(':')
        return f"{refseq_chromosome}:{g_dot_only}"
        
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
                    f.write(f"{v}\n")
                
                self._logger.info(f"Wrote {len(batch)} variant transcripts to {full_path}")
        

def _parse_args():
    parser = argparse.ArgumentParser(description='Read a list of variants and write them out in format suitable for submission to Mutalyzer''s Batch Processor')

    parser.add_argument('--variant_nomenclature', help='Variant nomenclature generated by workflow (csv)', required=True)
    parser.add_argument('--mutalyzer_cache', help="Local Mutalyzer cache (csv)", required=False)
    parser.add_argument('--out_dir', help='Directory to write batch files', required=True)

    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

def main():
    log_config = LogConfig().stdout_config
    log_config['loggers']['__main__']['level'] = 'INFO'
    logging.config.dictConfig(log_config)
    
    logging.debug("Hello")
    logging.info("Hi")
    
    args = _parse_args()
    
    f2b = FinalCsvToMutalyzerBatch(args.mutalyzer_cache)
    
    if args.mutalyzer_cache:
        logging.info("Using local Mutalyzer to filter variants we've already submitted")
    else:
        logging.info("Local Mutalyzer cache not provided for additional filtering")
    
    # Read and filter variants
    variants = f2b.get_variants(args.variant_nomenclature)
        
    # Write variants to batch files to be submitted to Mutalyzer's Batch Processor 
    f2b.write_transcripts(args.out_dir, variants)
    
if __name__ == '__main__':
    main()