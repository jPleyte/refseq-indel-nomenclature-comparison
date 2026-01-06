'''
Created on Jan 2, 2026

Generates a list of variants that are just downstream from insertion or deletion alignment gaps. 
   
@author: pleyte
'''
import logging.config
import argparse

import re
import pysam
from rinc.util import chromosome_map
import csv
from rinc.util.uta_db import UtaDb
from rinc.util.log_config import LogConfig

class FindGapVariants(object):
    '''
    Query the UTA database to find insertion or deletion alignment differences between hg19 and RefSeq transcript sequences 
    '''

    def __init__(self, uta_schema):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._uta_schema = uta_schema
    
    def _get_transcript_gaps(self):
        """
        Query the UTA database for hg19 transcripts that have insertion or deletion differences. 
        Cigar strings are limited to those that are pretty simple like 160=1D28= or 641=1I425=. 
        """
        cigar_regex = '^[0-9]+=[0-9]+[I|D][0-9]+=$'
        sql = f"""
                SELECT *
                  FROM {self._uta_schema}.tx_exon_aln_v 
                 WHERE alt_aln_method = 'splign'
                   AND tx_ac LIKE 'NM_%'
                   AND ord > 0
                   AND alt_ac IN (
                                  'NC_000001.10', 'NC_000002.11', 'NC_000003.11', 'NC_000004.11', 
                                  'NC_000005.9',  'NC_000006.11', 'NC_000007.13', 'NC_000008.10', 
                                  'NC_000009.11', 'NC_000010.10', 'NC_000011.9',  'NC_000012.11', 
                                  'NC_000013.10', 'NC_000014.8',  'NC_000015.9',  'NC_000016.9', 
                                  'NC_000017.10', 'NC_000018.9',  'NC_000019.9',  'NC_000020.10', 
                                  'NC_000021.8',  'NC_000022.10', 'NC_000023.10', 'NC_000024.9'
                              ) 
                   AND cigar ~ '{cigar_regex}'
                 LIMIT 10;
                """
        with UtaDb() as uta_db:
            transcript_gaps = uta_db.query(sql)
            self._logger.debug(f"Found {len(transcript_gaps)} transcript gaps")
            return transcript_gaps
        
    def _get_gap_coordinate(self, cigar: str, exon_start):
        """
        Use the exon start position and the cigar string to to find the genomic coordinate of the insertion or deletion  
        """
        # Split CIGAR (eg 641=1I425=) into components like [('641', '='), ('2', 'I'), ...]
        ops = re.findall(r'(\d+)([ID=])', cigar)
        current_genomic_location = exon_start
        
        for length, op in ops:
            if op == '=':
                current_genomic_location += int(length)
            elif op in 'ID':
                # We found the first Indel
                return current_genomic_location, op
            else:
                raise ValueError("Should not happen a")
        
        raise ValueError("Should not happen b")
                
    def find_gaps(self):
        """
        Query the UTA database for cigar strings with a single insertion or deletion.
        Get the genomic location of the insertion delete.
        Go 5 bases downstream, and if we're still on the exon keep the coordinate. 
        Return lis of these gaps.
        """
        gaps = []
        for t in self._get_transcript_gaps():
            exon_start = t['alt_start_i']
            exon_end = t['alt_end_i']
            strand = t['alt_strand']
            
            gap_coord, gap_type = self._get_gap_coordinate(t['cigar'], exon_start)
            
            # Figure out the coordinate that is five bases downstream
            step = 5 if strand == 1 else -5
            five_bases_downstream = gap_coord + step
            
            # Make sure the downstream position is still on the exon because we aren't interested in intronic changes
            if exon_start <= five_bases_downstream <= exon_end:
                t['five_bases_downstream'] = five_bases_downstream
                t['gap_type'] = 'deletion' if gap_type == 'D' else 'insertion'  
                gaps.append(t)

        return gaps

    def _get_variant_transcript(self, gap, chromosome, position, reference, alt):
        """
        Create a VariantTranscript object from the query result row 
        """
        variant_transcript = {
                              'chromosome': chromosome,
                              'position': position,
                              'reference': reference,
                              'alt': alt,
                              'cdna_transcript': gap['tx_ac'],
                              'strand': gap['alt_strand'],
                              'cigar': gap['cigar'],
                              'alt_start_i': gap['alt_start_i'],
                              'alt_end_i': gap['alt_end_i'],
                              'ord': gap['ord'],
                              'alt_ac': gap['alt_ac'],
                              'symbol': gap['symbol']                              
                              }
        
        return variant_transcript
        
    def get_variants(self, gaps: list[dict], fasta_file):
        """
        """
        fasta = pysam.FastaFile(fasta_file)
        variant_transcripts = []
        
        for row in gaps:
            chromosome = chromosome_map.get_ncbi(row['alt_ac'])
            
            # subtract one because pysam uses zero-based positions
            reference_base = fasta.fetch(reference=chromosome, start=row['five_bases_downstream']-1, end=row['five_bases_downstream'])
            
            arbitrary_different_base = { 'A': 'G', 'T': 'C', 'G': 'T', 'C': 'A' }
            variant_transcripts.append(self._get_variant_transcript(row, 
                                                                    chromosome, 
                                                                    row['five_bases_downstream'], 
                                                                    reference_base, 
                                                                    arbitrary_different_base[reference_base]))

        return variant_transcripts
    
    def write(self, out_filename, variants):
        """
        Write variant and gap information to csv file
        """
        headers = variants[0].keys()
        with open(out_filename, 'w', newline='') as output:
            writer = csv.DictWriter(output, fieldnames=headers)
            writer.writeheader()
            writer.writerows(variants)
            
        self._logger.info(f"Wrote {len(variants)} variants to {out_filename}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Query UTA for alignment differences')
    parser.add_argument("--version", action="version", version="0.0.1")
    parser.add_argument("--fasta", help="Geneome reference (fasta)", required=True)
    parser.add_argument("--uta_schema", help="UTA db schema", required=True)
    parser.add_argument("--out", help="output file (csv)", required=True)
    args = parser.parse_args()
    return args
    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    fgv = FindGapVariants(args.uta_schema)
    gaps = fgv.find_gaps()
    variants = fgv.get_variants(gaps, args.fasta)
    fgv.write(args.out, variants)
    
if __name__ == '__main__':
    main()