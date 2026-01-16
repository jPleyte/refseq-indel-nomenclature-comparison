'''
This module has the MytalyzerCache class which is a database of results 
received from Mutalyzer's Batch Processor (https://mutalyzer.nl/batchprocessor).
The main method is used to import batch files into the database file which is a csv
maintained in the git repository.
   
Created on Jan 11, 2026

@author: pleyte
'''
import argparse
import csv
import logging.config
import os
import re
import pandas as pd
from rinc.util.log_config import LogConfig

class MytalyzerCache(object):
    '''
    A local database of mutalyzer results
    '''
    def __init__(self, cache_db_file):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._cache_db_file = cache_db_file
        self._cache_db = self._load()
    
    def find(self, refseq_chromosome, g_dot, refseq_cdna_transcript):
        """
        Return result matching parameters or None if not found
        """
        results = self._cache_db[(self._cache_db['refseq_chromosome'] == refseq_chromosome) & 
                                 (self._cache_db['g_dot'] == g_dot) & 
                                 (self._cache_db['cdna_transcript'] == refseq_cdna_transcript)]
        
        if results.empty:
            # Try looking it up without transcript  
            results = self._cache_db[(self._cache_db['refseq_chromosome'] == refseq_chromosome) & 
                                 (self._cache_db['g_dot'] == g_dot)]
        
        if results.shape[0] > 1: 
            raise ValueError(f"Multiple rows matching: refseq_chromosome={refseq_chromosome}, g_dot={g_dot}, transcript={refseq_cdna_transcript}")
        elif results.shape[0] == 1:
            return results.iloc[0]
        else:
            return results
        
    def _load(self) -> dict:
        """
        Open and read in all rows in the local database
        """
        variant_schema = {
            'input_description': str,  
            'normalized': str,
            'status': str,
            'refseq_chromosome': str,
            'cdna_transcript': str,
            'protein_transcript': str,
            'g_dot': str,
            'c_dot': str,
            'p_dot': str,            
            'no_nomenclature': bool
        }
        
        # Return empty db if db file doesn't exist yet
        if not os.path.isfile(self._cache_db_file):
            self._logger.info(f"Local database {self._cache_db_file} does not exist yet.")
            cols = variant_schema.keys() 
            df = pd.DataFrame(columns=cols)
            return df.astype(variant_schema)        
        
        df = pd.read_csv(self._cache_db_file, dtype=variant_schema)
            
        self._logger.info(f"Read {df.shape[0]} rows from {self._cache_db_file}")
        return df
        
    def save_db(self):
        """
        Write all values to file        
        """
        self._cache_db.to_csv(self._cache_db_file, index=False)
        self._logger.info(f"Wrote {self._cache_db.shape[0]} results to {self._cache_db_file}")
    
    def _remove_equivalent_results(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        We may end up with duplicate rows as a result of requesting the same thing different ways. 
        Example: The g. request "NC_000016.9:g.3602211C>T" and the c. request "NC_000016.9(NM_178844.4):c.2336G>A" produce the same result.
            When that happens we must drop one in order to maintain row uniqueness among the rest of the fields. 
        """
        mask = df['input_description'].str.contains('g.', na=False)
        filtered_df = df[mask]
        for index, row in filtered_df.iterrows():
            refseq_chromosome = row['refseq_chromosome']
            g_dot = row['g_dot']
            cdna_transcript = row['cdna_transcript']
            
            duplicates = df[(df['refseq_chromosome'] == refseq_chromosome) & 
                            (df['g_dot'] == g_dot) & 
                            (df['cdna_transcript'] == cdna_transcript) &
                            (df['input_description'].str.contains('c.'))
                            ]
            
            if not duplicates.empty:
                assert duplicates.shape[0] == 1, f"When there is a duplicate there should be only one: chr={refseq_chromosome}, g.={g_dot}, transcript={cdna_transcript}"
                
                self._logger.info(f"Dropping c. equivalent of of {row['input_description']}")
                df.drop(duplicates.index, inplace=True)
        
        return df
                
    def import_results(self, batch_result_file):
        """
        Read the csv containing mutalyzer results 
        """
        new_results = []
        with open(batch_result_file, mode='r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                mr = self._read_mutalyzer_result(row)
                new_results.append(mr)
        
        new_df = pd.DataFrame(new_results)
        
        old_size = self._cache_db.shape[0]
        
        
        self._cache_db = self._remove_equivalent_results(pd.concat([self._cache_db, new_df]).drop_duplicates(ignore_index=True))
        
                
        new_size = self._cache_db.shape[0]
        
        self._logger.info(f"Imported {new_size-old_size} new results. New size is {new_size}")
                
    def _read_mutalyzer_result(self, row):
        """
        Extract values from a row read from the Batch Processor csv
        """            
        mutalyzer_result = {}
        mutalyzer_result['input_description'] = row['Input description']
        mutalyzer_result['normalized'] = row['Normalized']
        mutalyzer_result['status'] = row['Status']

        if row['Status'] == "Failed":
            # When status=Faild it'probably s because the reference values were incorrect.
            # Should we set mutalyzer_result['no_nomenclature'] = True?
            mutalyzer_result['no_nomenclature'] = False 
            return mutalyzer_result
                
        # Get the g. which comes from different fields depending on the original request type
        mutalyzer_result['refseq_chromosome'], mutalyzer_result['g_dot'] = self.get_g_dot_parts(row)
        
        # If the request was a g. string then parse c. and cdna transcript from "DNA transcript" field. May return three None values
        # When the request was a c. then it sets DNA Transcript to "N/A" and we just parse our c. out of the Normalized field. 
        dna_genomic = row['Normalized'] if ("c." in row['Input description'] and row['DNA transcript'] == "N/A") else row['DNA transcript']         
        refseq_chromosome_d, cdna_transcript, c_dot = self._get_dna_transcript_parts(dna_genomic)
        if not (refseq_chromosome_d is None and cdna_transcript is None and c_dot is None):
            assert mutalyzer_result['refseq_chromosome'] == refseq_chromosome_d, "Normalized chromosome and DNA chromosome should be the same"
            mutalyzer_result['cdna_transcript'] = cdna_transcript
            mutalyzer_result['c_dot'] = c_dot
        
        # Parse p. and protein transcript from "Protein" field. May return three None values
        refseq_chromosome_p, protein_transcript, p_dot = self._get_protein_parts(row['Protein'])
        if not (refseq_chromosome_p is None and protein_transcript is None and p_dot is None):
            assert refseq_chromosome_p == refseq_chromosome_d, f"Normalized chromosome and protein chromosome should be the same: {refseq_chromosome_p} != {refseq_chromosome_d}. See {row['Input description']}"
            mutalyzer_result['protein_transcript'] = protein_transcript
            mutalyzer_result['p_dot'] = self._get_pdot_without_parenthesis(p_dot)
        
        if refseq_chromosome_d is None and cdna_transcript is None and c_dot is None and refseq_chromosome_p is None and protein_transcript is None and p_dot is None:                    
            mutalyzer_result['no_nomenclature'] = True
        else: 
            mutalyzer_result['no_nomenclature'] = False
        
        return mutalyzer_result
    
    def get_g_dot_parts(self, row: dict) -> (str, str):
        """
        Parse "g." and refseq chromosome (NC_) from the "DNA genomic" or "Normalized" field
        """
        if row['DNA genomic'] != 'N/A':
            return self._get_normalized_parts(row['DNA genomic'])
        else:
            return self._get_normalized_parts(row['Normalized'])
                      
    def _get_pdot_without_parenthesis(self, p_dot: str) -> str:
        """
        Remove parenthesis from p. (even though it is wrong to do so)
        """
        match = re.match(r"p\.\((.*)\)", p_dot)
        if match:
            return f"p.{match.group(1)}"
        return p_dot
        
    def _get_normalized_parts(self, normalized: str):
        """
        Read the Batch Processor's "Normalized" field and return the chromosome and g.
        normalized string should look like "NC_000012.11:g.9994422G>A" 
        """
        assert normalized, "Normalized value must not be empty"
        ref_chr, g_dot = normalized.split(':')
        
        assert ref_chr.startswith("NC_"), "Refseq chromosome must start with NC"
        assert g_dot.startswith("g."), "g. must start with g."
        
        return ref_chr, g_dot

    def _get_dna_transcript_parts(self, dna_transcript: str):
        """
        Read the Batch Processor's "DNA transcript" field and return the chormosome, transcript, and g.
        DNA transcript string should look like "NC_000016.9(NM_178844.4):c.2336G>A" or "N/A"
        """
        if dna_transcript == 'N/A':
            return None, None, None
    
        pattern = r"(?P<chrom>[^(]+)\((?P<transcript>[^)]+)\):(?P<c_dot>.+)"
        match = re.search(pattern, dna_transcript)
        
        if match:            
            chromosome = match.group('chrom')
            assert chromosome.startswith("NC_"), "Refseq chromosome must start with NC"
            
            transcript = match.group('transcript')
            assert transcript.startswith("NM_"), "cDNA transcript must start with NM"
            
            c_dot = match.group('c_dot')
            assert c_dot.startswith("c."), "c. must start with c."
            
            return chromosome, transcript, c_dot
        else:
            raise ValueError(f"DNA Transcript field does not match expected pattern: {dna_transcript}") 
    
    def _get_protein_parts(self, protein: str):
        """
        Read the Batch Processor's "Protein" field and return the chromosome, protein transcript, and p.
        Protein string should look like "NC_000016.9(NP_849172.2):p.(Ser779Asn)" or "N/A"
        """
        if protein == 'N/A':
            return None, None, None
        
        pattern = r"(?P<chrom>[^(]+)\((?P<transcript>[^)]+)\):(?P<p_dot>p\..+)"
        match = re.search(pattern, protein)
        
        if match:
            chromosome = match.group('chrom')
            assert chromosome.startswith("NC_"), "Refseq chromosome must start with NC"
            
            transcript = match.group('transcript')
            assert transcript.startswith("NP_"), "Protein transcript must start with NP"
            
            p_dot = match.group('p_dot')
            assert p_dot.startswith("p."), "p. must start with p."
            
            return chromosome, transcript, self._get_pdot_without_parenthesis(p_dot)            
        else:
            raise ValueError(f"Protein field does not match expected pattern: {protein}")


def _parse_args():
    parser = argparse.ArgumentParser(description='Import mutalyzer Batch Processor results into the local cache')
    parser.add_argument("--version", action="version", version="0.0.1")    
    parser.add_argument("--batch_file", help="mutalyzer Batch Processor results to import (csv)", required=True)
    parser.add_argument("--mutalyzer_cache", help="Local mutalyzer cache where (csv)", required=True)
    
    args = parser.parse_args()
    return args    

def main():
    logging.config.dictConfig(LogConfig().stdout_config)

    args = _parse_args()
    
    mc = MytalyzerCache(args.mutalyzer_cache)
    mc.import_results(args.batch_file)
    mc.save_db()
    
if __name__ == '__main__':
    main()        
        