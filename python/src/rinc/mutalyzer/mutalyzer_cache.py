'''
This module has the MytalyzerCache class which is a database of results 
received from Mutalyzer's Batch Processor (https://mutalyzer.nl/batchprocessor).
The main method is used to import batch files into the database file which is a csv
maintained in the git repository.
   
Created on Jan 11, 2026

@author: pleyte
'''
import argparse
from collections import Counter
import csv
from dataclasses import dataclass, asdict
import logging.config
import os
import re

from rinc.util.log_config import LogConfig


@dataclass(slots=True)
class MutilzizerResult():
    input_description: str
    
    normalized: str = None
    refseq_chromosome: str = None
    cdna_transcript: str = None
    protein_transcript: str = None
    g_dot: str = None
    c_dot: str = None
    p_dot: str = None
    
    # When true, indicates that this variant has been queried and normalied but no transcripts found at this location
    is_empty: bool = False
    
    def __post_init__(self):
        # Fix the is_empty_ boolean type if it was read as a string
        if isinstance(self.is_empty, str):
            self.is_empty = (self.is_empty.lower() == "true")

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
    
    def find(self, refseq_chromosome, g_dot, refseq_cdna_transcript) -> MutilzizerResult:
        """
        Return result matching parameters or None if not found
        """
        return self._cache_db.get(self._get_db_key(refseq_chromosome, g_dot, refseq_cdna_transcript))
        
    def _load(self) -> dict:
        """
        Open and read in all rows in the local database
        """
        # Return empty db if db file doesn't exist yet
        if not os.path.isfile(self._cache_db_file):
            self._logger.info(f"Local database {self._cache_db_file} does not exist yet.")
            return {}
        
        with open(self._cache_db_file, 'r', encoding='utf-8') as f:
            db = {}
            reader = csv.DictReader(f)
            for row in reader:
                mr =  MutilzizerResult(**row)
                key = self._get_db_key_mr(mr)
                assert key not in db
                db[key] = mr
            
            self._logger.info(f"Read {len(db)} from {self._cache_db_file}")
            return db
        
    def save_db(self):
        """
        Write all values to file        
        """
        with open(self._cache_db_file, 'w', newline='', encoding='utf-8') as f:
            # Use the class fields as header names
            writer = csv.DictWriter(f, fieldnames=MutilzizerResult.__annotations__.keys())
            writer.writeheader()
            for x in self._cache_db.values():
                writer.writerow(asdict(x))
        
        self._logger.info(f"Wrote {len(self._cache_db)} results to {self._cache_db_file}")
        
    def _get_db_key_mr(self, mr: MutilzizerResult):
        """
        Return the database key for MutilzizerResult 
        """
        if mr.is_empty:
            return self._get_db_key(mr.refseq_chromosome, mr.g_dot, "noTranscript")
        else:
            return self._get_db_key(mr.refseq_chromosome, mr.g_dot, mr.cdna_transcript)
    
    def _get_db_key(self, refseq_chromosome, g_dot, refseq_cdna_transcript):
        """
        Return the database key constructed from the given fields
        """
        transcript = "noTranscript" if refseq_cdna_transcript is None else refseq_cdna_transcript  
        assert refseq_chromosome.startswith("NC")
        assert g_dot.startswith("g."), f"Expecting g. that starts with g. but got {g_dot}"
        assert transcript == 'noTranscript' or refseq_cdna_transcript.startswith("NM") 
                
        return f"{refseq_chromosome}-{g_dot}-{transcript}"
    
    def _are_equal(self, a:MutilzizerResult, b: MutilzizerResult):
        if not a or not b:
            raise ValueError(f"left or right side of comparison is missing: {a} == {b}")
        if a.input_description != b.input_description:
            return False
        
        if a.normalized != b.normalized:
            return False
        
        if a.is_empty and b.is_empty:
            return True
        
        return (a.g_dot == b.g_dot and
                a.refseq_chromosome == b.refseq_chromosome and
                a.c_dot == b.c_dot and
                a.cdna_transcript == b.cdna_transcript and
                a.p_dot == b.p_dot and
                a.protein_transcript == b.protein_transcript)
                 
        
    def import_results(self, batch_result_file):
        """
        Read the csv containing mutalyzer results 
        """
        count = Counter()
        with open(batch_result_file, mode='r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                if row['Status'] == "Failed":
                    print(f"Encountered failed variant: {row['Input description']}")
                    count['failed'] += 1
                    continue
                
                mr = self._read_row(row)

                db_key = self._get_db_key_mr(mr)
                if db_key in self._cache_db:
                    if self._are_equal(mr, self._cache_db.get(db_key)):
                        count['import_already_exists'] += 1
                        self._logger.info(f"Record already exists, skipping: {mr.input_description}")
                    else:
                        raise ValueError(f"Local values and incoming values are different for: {mr.input_description}")
                else:
                    self._cache_db[db_key] = mr
                    count['import_new_record'] += 1
        
        self._logger.info(f"Imported {batch_result_file}: {count}")
                
    def _read_row(self, row):
        """
        Extract values from a row read from the Batch Processor csv
        """
        mr = MutilzizerResult(row['Input description'])
        mr.normalized = row['Normalized']
                
        # Parse "g." and refseq chromosome (NC_) from Normalized field
        refseq_chromosome_n, g_dot = self._get_normalized_parts(row['Normalized'])
        mr.refseq_chromosome = refseq_chromosome_n
        mr.g_dot = g_dot
        
        # Parse c. and cdna transcript from "DNA transcript" field. May return three None values
        refseq_chromosome_d, cdna_transcript, c_dot = self._get_dna_transcript_parts(row['DNA transcript'])
        if not (refseq_chromosome_d is None and cdna_transcript is None and c_dot is None):
            assert refseq_chromosome_n == refseq_chromosome_d, "Normalized chromosome and DNA chromosome should be the same"
            mr.cdna_transcript = cdna_transcript
            mr.c_dot = c_dot
        
        # Parse p. and protein transcript from "Protein" field. May return three None values
        refseq_chromosome_p, protein_transcript, p_dot = self._get_protein_parts(row['Protein'])
        if not (refseq_chromosome_p is None and protein_transcript is None and p_dot is None):
            assert refseq_chromosome_p == refseq_chromosome_d, "Normalized chromosome and protein chromosome should be the same"
            mr.protein_transcript = protein_transcript
            mr.p_dot = self._get_pdot_without_parenthesis(p_dot)
        
        if refseq_chromosome_d is None and cdna_transcript is None and c_dot is None and refseq_chromosome_p is None and protein_transcript is None and p_dot is None:                    
            mr.is_empty = True
        
        return mr
    
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
        try:
            ref_chr, g_dot = normalized.split(':')
        except ValueError as e:
            print(e)
        
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
        