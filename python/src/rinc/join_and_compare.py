'''
Created on Jan 5, 2026
Reader Gather gap, hgvs, and annovar information; perform comparison, and write results to csv.
  
@author: pleyte
'''
import argparse
import csv
import logging.config
from rinc.util.log_config import LogConfig
from rinc.variant_transcript import VariantTranscript
import pandas as pd
import pandas
from functools import reduce

class JoinAndCompare(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        
    
    def get_gap_variants(self, gap_variants_file) -> pd.DataFrame:
        """
        Return gap information keyed by variant and transcript 
        """
        return pd.read_csv(gap_variants_file) 
        
    
    def get_hgvs_nomenclature(self, hgvs_nomenclature_file) -> pd.DataFrame:
        """
        Return hgvs nomenclature in a dict keyed by variant and transcript
        """
        return pd.read_csv(hgvs_nomenclature_file)

    def get_annovar_nomenclature(self, annovar_nomenclature_file) -> pd.DataFrame:
        """
        Return annovar nomenclature in a dict keyed by variant and transcript
        """
        return pd.read_csv(annovar_nomenclature_file)

    def get_comparison(self, gaps_df: pd.DataFrame, hgvs_df: pd.DataFrame, annovar_df: pd.DataFrame):
        """
        """
        create_mask = lambda df, c, p, r, a, t: (
            (df['chromosome'] == c) & 
            (df['position'] == p) & 
            (df['reference'] == r) & 
            (df['alt'] == a) &
            (df['cdna_transcript'] == t))
        
        comparisons = []
        
        for row in gaps_df[['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']].itertuples(index=False):
            chromosome = row.chromosome
            position = row.position
            reference = row.reference
            alt = row.alt
            cdna_transcript = row.cdna_transcript
            
            hgvs_row = hgvs_df[create_mask(hgvs_df, chromosome, position, reference, alt, cdna_transcript)]
            annovar_row = annovar_df[create_mask(annovar_df, chromosome, position, reference, alt, cdna_transcript)]
            
            comparisons.append(self._get_comparison(chromosome, position, reference, alt, cdna_transcript, hgvs_row, annovar_row))
        
        return pd.DataFrame(comparisons)

    def _get_comparison(self, chromosome, position, reference, alt, cdna_transcript, hgvs_row, annovar_row):
        """
        """
        comparison = {'chromosome': chromosome,
                       'position': position,
                       'reference': reference,
                       'alt': alt,
                       'cdna_transcript': cdna_transcript}
        
        if hgvs_row.empty: 
            raise ValueError("hgvs should always have a result")
        elif annovar_row.empty:
            comparison['differences'] = 'hgvs_only'
        else:            
            c_dot_are_different = hgvs_row['hu.c_dot'].item() != annovar_row['annovar.c_dot'].item()
            p_dot_are_different = hgvs_row['hu.p_dot1'].item() != annovar_row['annovar.p_dot1'].item()
            exon_are_different = hgvs_row['hu.exon'].item() != annovar_row['annovar.exon'].item()
            difference_count = c_dot_are_different + p_dot_are_different + exon_are_different
            comparison['differences'] = difference_count
            
        return comparison
    
    def write(self, out_file, gaps_df: pd.DataFrame, hgvs_df: pd.DataFrame, annovar_df: pd.DataFrame, comparison_df: pd.DataFrame):
        """
        """
        dataframes = [gaps_df, hgvs_df, annovar_df, comparison_df]
        join_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        
        # Perform outter join on variant+transcript fields 
        merged_df = reduce(lambda left, right: pd.merge(left, right, on=join_cols, how='outer'), dataframes)
        merged_df.to_csv(out_file, index=False, encoding='utf-8')
        self._logger.info(f"Wrote data frame with {merged_df.shape[0]} rows and {merged_df.shape[1]} columns  to {out_file}")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Read Annovar multianno file, extract values and write to new csv')

    parser.add_argument('--gap_variants', help='Gap information and variants (csv)', required=True)
    parser.add_argument('--hgvs_nomenclature', help='hgvs/uta values (csv)', required=True)
    parser.add_argument('--annovar_nomenclature', help='Annovar values (csv)', required=True)    
    

    parser.add_argument("--out", help="output file (csv)", required=True)

    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    jc = JoinAndCompare()
    
    gaps = jc.get_gap_variants(args.gap_variants)
    
    hgvs = jc.get_hgvs_nomenclature(args.hgvs_nomenclature)
    
    annovar = jc.get_annovar_nomenclature(args.annovar_nomenclature)    

    comparison = jc.get_comparison(gaps, hgvs, annovar)
    
    jc.write(args.out, gaps, hgvs, annovar, comparison)
    

if __name__ == '__main__':
    main()
        