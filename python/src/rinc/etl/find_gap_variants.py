'''
Created on Jan 2, 2026

Generates a list of variants that are just downstream from insertion or deletion alignment gaps. 
   
@author: pleyte
'''
import logging.config
import argparse

import re
from rinc.util import chromosome_map
import csv
from rinc.util.uta_db import UtaDb
from rinc.util.log_config import LogConfig
from collections import Counter
from rinc.util.tx_eff_pysam import PysamTxEff
from rinc.util.vcf_to_gdot import get_gdot
from rinc.variant_transcript import VariantTranscript

def get_variants(variants_file: str):
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
    return variants

class FindGapVariants(object):
    '''
    Query the UTA database to find insertion or deletion alignment differences between hg19 and RefSeq transcript sequences 
    '''

    def __init__(self, uta_schema, fasta_file):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._uta_schema = uta_schema
        self._pysam_tx = PysamTxEff(fasta_file)
    
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
                  LIMIT 500;
                """
                          
        with UtaDb() as uta_db:
            transcript_gaps = uta_db.query(sql)
            self._logger.debug(f"Found {len(transcript_gaps)} transcript gaps")
            return transcript_gaps
        
    def _get_gap_coordinate(self, strand: int, cigar: str, exon_start: int, exon_end: int):
        """
        Use the exon start position and the cigar string to to find the genomic coordinate of the insertion or deletion  
        """
        # Split CIGAR (eg 641=1I425=) into components like [('641', '='), ('2', 'I'), ...]
        ops = re.findall(r'(\d+)([ID=])', cigar)
        
        assert strand == -1 or strand == 1, f"Strand must be -1 or 1 {strand}"
        assert exon_start < exon_end, f"Exon start is after exon end: {exon_start}, {exon_end}"
        
        # Transcripts on positive strand start at 
        genomic_position = exon_start if strand == 1 else exon_end
        
        for length, op in ops:
            if op == '=':
                genomic_position += int(length) * strand
            elif op in 'ID':
                # We found the first Indel
                return genomic_position, op
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
        find_gap_counter = Counter()
        
        gaps = []
        for t in self._get_transcript_gaps():
            exon_start = t['alt_start_i']
            exon_end = t['alt_end_i']
            strand = t['alt_strand']
            
            gap_coord, gap_type = self._get_gap_coordinate(int(strand), t['cigar'], exon_start, exon_end)
            
            # jdebug: do we also need to deal with deletions or insertions somehow? 
            
            # Figure out the coordinate that is five bases downstream
            step = 5 if strand == 1 else -5
            five_bases_downstream = gap_coord + step
            
            # Make sure the downstream position is still on the exon because we aren't interested in intronic changes
            if exon_start <= five_bases_downstream <= exon_end:
                t['five_bases_downstream'] = five_bases_downstream
                t['gap_type'] = 'deletion' if gap_type == 'D' else 'insertion'  
                gaps.append(t)
                find_gap_counter['downstream_position_found'] += 1
            else: 
                find_gap_counter['downstream_beyond_exon'] += 1

        self._logger.info(f"Downstream position results: {find_gap_counter}")
        return gaps

    def _get_row_details(self, gap, chromosome, position, reference, alt):
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
                              'symbol': gap['symbol'],      
                              'g_dot': self._get_g_dot(chromosome, position, reference, alt)
                              }
        
        return variant_transcript
        
    def _get_g_dot(self, chromosome, position, reference, alt):
        """
        Return the variant using g. nomenclature 
        """
        gdot = get_gdot(chromosome, position, reference, alt, self._pysam_tx)
        refseq_chromosome = chromosome_map.get_refseq(chromosome)
        return refseq_chromosome + ':g.' + gdot
    
        
    def get_variants(self, gaps: list[dict]):
        """
        """
        variant_transcripts = []
        
        for row in gaps:
            chromosome = chromosome_map.get_ncbi(row['alt_ac'])
            
            # subtract one because pysam uses zero-based positions
            reference_base = self._pysam_tx.direct_query(chromosome, row['five_bases_downstream']-1, row['five_bases_downstream'])
            
            arbitrary_different_base = { 'A': 'G', 'T': 'C', 'G': 'T', 'C': 'A' }
            
            if reference_base == 'N':
                self._logger.debug(f"Invalid variant having base 'N': {chromosome} {row['alt_ac']} {row['five_bases_downstream']}")
                continue
            variant_transcripts.append(self._get_row_details(row, 
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
    
    fgv = FindGapVariants(args.uta_schema, args.fasta)
    gaps = fgv.find_gaps()
    variants = fgv.get_variants(gaps)
    fgv.write(args.out, variants)
    
if __name__ == '__main__':
    main()