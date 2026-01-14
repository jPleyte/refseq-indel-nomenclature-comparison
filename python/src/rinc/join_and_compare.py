'''
Created on Jan 5, 2026
Gather gap, hgvs, and annovar information; perform comparison, and write results to csv.
  
@author: pleyte
'''
import argparse
import logging.config
from rinc.util.log_config import LogConfig
import pandas as pd
from functools import reduce
from itertools import combinations

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
        return pd.read_csv(gap_variants_file, dtype=str) 
        
    
    def get_hgvs_nomenclature(self, hgvs_nomenclature_file) -> pd.DataFrame:
        """
        Return the hgvs nomenclature dataframe
        """
        hgvs_df = pd.read_csv(hgvs_nomenclature_file, dtype=str)
        return hgvs_df

    def get_annovar_nomenclature(self, annovar_nomenclature_file) -> pd.DataFrame:
        """
        Return annovar nomenclature dataframe
        """
        annovar_df = pd.read_csv(annovar_nomenclature_file, dtype=str)
        return annovar_df 

    def get_snpeff_nomenclature(self, snpeff_nomenclature_file) -> pd.DataFrame:
        """
        Return the SnpEff nomenclature dataframe
        """
        annovar_df = pd.read_csv(snpeff_nomenclature_file, dtype=str)
        return annovar_df
    
    def get_mutalyzer_nomenclature(self, mutalyzer_nomenclature_file) -> pd.DataFrame:
        """
        Return the mutalyzer nomenclature dataframe
        """
        mutalyzer_df = pd.read_csv(mutalyzer_nomenclature_file, dtype=str)
        return mutalyzer_df
    
    def get_vep_nomenclature(self, vep_nomenclature_file) -> pd.DataFrame:
        """
        Return the mutalyzer nomenclature dataframe
        """
        vep_nomenclature_df = pd.read_csv(vep_nomenclature_file, dtype=str)
        return vep_nomenclature_df


    def get_comparison_df(self, gaps_df: pd.DataFrame, 
                          hgvs_df: pd.DataFrame, 
                          annovar_df: pd.DataFrame, 
                          snpeff_df: pd.DataFrame, 
                          mutalyzer_df: pd.DataFrame, 
                          vep_refseq_nomenclature_df: pd.DataFrame, 
                          vep_hg19_nomenclature_df: pd.DataFrame):
        """
        Compare dataframes with each other and return a dataframe with comparison score         
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
            
            # Ideally the mask would select exactly one row but i've seen duplicates come from annovar 
            # when the same transcript is found on multiple genes, so it is necessary to use ``.drop_duplicates()``
            annovar_row = annovar_df[create_mask(annovar_df, chromosome, position, reference, alt, cdna_transcript)].drop_duplicates()
            
            snpeff_row = snpeff_df[create_mask(snpeff_df, chromosome, position, reference, alt, cdna_transcript)]

            mutalyzer_row = mutalyzer_df[create_mask(mutalyzer_df, chromosome, position, reference, alt, cdna_transcript)]
            
            vep_refseq_row = vep_refseq_nomenclature_df[create_mask(vep_refseq_nomenclature_df, chromosome, position, reference, alt, cdna_transcript)]
            
            vep_hg19_row = vep_hg19_nomenclature_df[create_mask(vep_hg19_nomenclature_df, chromosome, position, reference, alt, cdna_transcript)]
            
            if hgvs_row.shape[0] > 1:
                raise ValueError(f"hgvs nomenclature has multiple rows for: {chromosome, position, reference, alt, cdna_transcript}")
            elif annovar_row.shape[0] > 1: 
                raise ValueError(f"annovar nomenclature has multiple rows for: {chromosome, position, reference, alt, cdna_transcript}")
            elif snpeff_row.shape[0] > 1:
                raise ValueError(f"snpeff nomenclature has multiple rows for: {chromosome, position, reference, alt, cdna_transcript}")
            elif mutalyzer_row.shape[0] > 1:
                raise ValueError(f"mutalyzer nomenclature has multiple rows for: {chromosome, position, reference, alt, cdna_transcript}")
            elif vep_refseq_row.shape[0] > 1:
                raise ValueError(f"vep_refseq nomenclature has multiple rows for: {chromosome, position, reference, alt, cdna_transcript}")
            elif vep_hg19_row.shape[0] > 1:
                raise ValueError(f"vep_hg19 nomenclature has multiple rows for: {chromosome, position, reference, alt, cdna_transcript}")
            
            comparisons.append(self._get_comparison(chromosome, position, reference, alt, cdna_transcript, 
                                                    hgvs_row, annovar_row, snpeff_row, mutalyzer_row,
                                                    vep_refseq_row, vep_hg19_row))
        
        return pd.DataFrame(comparisons).astype({'pairwise_score': float})

    def _get_comparison(self, chromosome, position, reference, alt, cdna_transcript, 
                        hgvs_row: pd.DataFrame, 
                        annovar_row: pd.DataFrame, 
                        snpeff_row: pd.DataFrame, 
                        mutalyzer_row: pd.DataFrame,
                        vep_refseq_row: pd.DataFrame, 
                        vep_hg19_row: pd.DataFrame) -> float:
        """
        Calculate pairwise comparison score for c., p., and exon from all sources.
        """
        comparison = {'chromosome': chromosome,
                       'position': position,
                       'reference': reference,
                       'alt': alt,
                       'cdna_transcript': cdna_transcript}

        if hgvs_row.empty: 
            raise ValueError(f"hgvs should always have a result: {comparison}")
        
        # c. values for comparison
        c_dot_values = [ hgvs_row['hu.c_dot'].item() ]
        if not annovar_row.empty:
            c_dot_values.append(annovar_row['annovar.c_dot'].item())
        if not snpeff_row.empty:
            c_dot_values.append(snpeff_row['snpeff.c_dot'].item())
        if not mutalyzer_row.empty:
            c_dot_values.append(mutalyzer_row['mut.c_dot'].item())
        if not vep_refseq_row.empty:
            c_dot_values.append(vep_refseq_row['vep.refseq.c_dot'].item())
        if not vep_hg19_row.empty:
            c_dot_values.append(vep_hg19_row['vep.hg19.c_dot'].item())
             
        # p. values for comparison
        p_dot_values = [ hgvs_row['hu.p_dot1'].item() ]
        if not annovar_row.empty:
            p_dot_values.append(annovar_row['annovar.p_dot1'].item())
        if not snpeff_row.empty:
            p_dot_values.append(snpeff_row['snpeff.p_dot1'].item())
        if not mutalyzer_row.empty:
            p_dot_values.append(mutalyzer_row['mut.p_dot1'].item())
        if not vep_refseq_row.empty:
            p_dot_values.append(vep_refseq_row['vep.refseq.p_dot1'].item())
        if not vep_hg19_row.empty:
            p_dot_values.append(vep_hg19_row['vep.hg19.p_dot1'].item())
                     
        # exon vlaues for comparison
        exon_values = [ hgvs_row['hu.exon'].item() ]
        if not annovar_row.empty:
            exon_values.append(annovar_row['annovar.exon'].item())
        if not snpeff_row.empty:
            exon_values.append(snpeff_row['snpeff.exon'].item())
        if not vep_refseq_row.empty:
            exon_values.append(vep_refseq_row['vep.refseq.exon'].item())
        if not vep_hg19_row.empty:
            exon_values.append(vep_hg19_row['vep.hg19.exon'].item())
        
        c_pairs  = list(combinations(c_dot_values, 2))
        p_pairs = list(combinations(p_dot_values, 2))
        e_pairs = list(combinations(exon_values, 2))
        
        c_matches = sum(1 for a, b in c_pairs if a == b)
        p_matches = sum(1 for a, b in p_pairs if a == b)
        e_matches = sum(1 for a, b in e_pairs if a == b)
        
        c_score = -1 if len(c_pairs) == 0 else (c_matches / len(c_pairs))
        p_score = -1 if len(p_pairs) == 0 else (p_matches / len(p_pairs))
        e_score = -1 if len(e_pairs) == 0 else (e_matches / len(e_pairs))
        
        if c_score == -1 and p_score == -1 and e_score == -1:
            # a score of -3 means hgvs/uta is the only source providing a value 
            final_score = -3
        else:
            score_sum = 0
            divisor = 0
            
            if c_score >= 0:
                score_sum += c_score
                divisor += 1
            
            if p_score >= 0:
                score_sum += p_score
                divisor += 1
            
            if e_score >= 0:
                score_sum += e_score
                divisor += 1            
        
            final_score = (c_score + p_score + e_score) / divisor
        
        comparison['pairwise_score'] = round(final_score, 3)
            
        return comparison
    
    def write(self, 
              out_file, 
              gaps_df: pd.DataFrame, 
              hgvs_df: pd.DataFrame, 
              annovar_df: pd.DataFrame, 
              snpeff_df: pd.DataFrame, 
              mutalyzer_df: pd.DataFrame,
              vep_refseq_nomenclature_df: pd.DataFrame, 
              vep_hg19_nomenclature_df: pd.DataFrame,
              comparison_df: pd.DataFrame):
        """
        Join all the dataframs on the variant+transcript fields and write out a csv with all fields.  
        """
        dataframes = [gaps_df, hgvs_df, annovar_df, snpeff_df, mutalyzer_df, vep_refseq_nomenclature_df, vep_hg19_nomenclature_df, comparison_df]
        join_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        
        # Perform outter join on variant+transcript fields 
        merged_df = reduce(lambda left, right: pd.merge(left, right, on=join_cols, how='outer'), dataframes)
        
        # There are a lot of fields in the dataframe. Move the important ones up front 
        front_columns = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript', 
                         'hu.exon', 'hu.tfx_exon', 'annovar.exon', 'snpeff.exon', 'vep.refseq.exon', 'vep.hg19.exon',
                         'hu.c_dot', 'annovar.c_dot', 'snpeff.c_dot', 'mut.c_dot', 'vep.refseq.c_dot', 'vep.hg19.c_dot',
                         'hu.p_dot1', 'annovar.p_dot1', 'snpeff.p_dot1', 'mut.p_dot1', 'vep.refseq.p_dot1', 'vep.hg19.p_dot1',
                         'pairwise_score', 
                         'g_dot', 'hu.g_dot', 'mut.g_dot', 'vep.refseq.g_dot', 'vep.hg19.g_dot', 
                         'hu.p_dot3', 'snpeff.p_dot3', 'mut.p_dot3', 'vep.refseq.p_dot3', 'vep.hg19.p_dot3',
                         'hu.protein_transcript', 'mut.protein_transcript', 'vep.refseqprotein_transcript', 'vep.hg19protein_transcript',
                         'symbol', 'hu.gene', 'annovar.gene', 'snpeff.gene', 'vep.refseq.gene', 'vep.hg19.gene',
                         'strand', 'hu.strand', 'ord',
                         'annovar.genomic_region_type', 'snpeff.genomic_region_type', 'vep.refseq.grt', 'vep.hg19.grt',
                         'annovar.protein_variant_type', 'ann_effect_raw', 'vep.refseq.pvt', 'vep.hg19.pvt', 'ann_biotype_raw', 
                         'hu.note', 'annovar.note',
                         'vep.refseq.REFSEQ_MATCH', 'vep.hg19.REFSEQ_MATCH', 'vep.refseq.BAM_EDIT']
        
        # Gather a list of all columsn other than the front_columns
        other_columns = [c for c in merged_df.columns if c not in front_columns]
        
        ordered_df = merged_df[front_columns + other_columns]
        
        sorted_df = ordered_df.sort_values(by=['pairwise_score','annovar.p_dot1', 'chromosome', 'position', 'reference', 'alt', 'cdna_transcript'], 
                                           ascending=[False, True, True, True, True, True, True])
        
        sorted_df.to_csv(out_file, index=False, encoding='utf-8')
        
        self._logger.info(f"Wrote data frame with {merged_df.shape[0]} rows and {merged_df.shape[1]} columns  to {out_file}")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Read Annovar multianno file, extract values and write to new csv')

    parser.add_argument('--gap_variants', help='Gap information and variants (csv)', required=True)
    parser.add_argument('--hgvs_nomenclature', help='hgvs/uta values (csv)', required=True)
    parser.add_argument('--annovar_nomenclature', help='Annovar values (csv)', required=True)    
    parser.add_argument('--snpeff_nomenclature', help='SnpEff values (csv)', required=True)
    parser.add_argument('--mutalyzer_nomenclature', help='mutalyzer values (csv)', required=True)
    parser.add_argument('--vep_refseq_nomenclautre', help='mutalyzer values (csv)', required=True)
    parser.add_argument('--vep_hg19_nomenclature', help='mutalyzer values (csv)', required=True)

    parser.add_argument("--out", help="output file (csv)", required=True)

    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    jc = JoinAndCompare()
    
    gaps_df = jc.get_gap_variants(args.gap_variants)
    
    hgvs_df = jc.get_hgvs_nomenclature(args.hgvs_nomenclature)
    
    annovar_df = jc.get_annovar_nomenclature(args.annovar_nomenclature)
    
    snpeff_df = jc.get_snpeff_nomenclature(args.snpeff_nomenclature)
    
    mutalyzer_df = jc.get_mutalyzer_nomenclature(args.mutalyzer_nomenclature)
    
    vep_refseq_nomenclature_df = jc.get_vep_nomenclature(args.vep_refseq_nomenclautre)
    
    vep_hg19_nomenclature_df = jc.get_vep_nomenclature(args.vep_hg19_nomenclature)

    # Compare exon, c., and p. from all datasources. 
    comparison_df = jc.get_comparison_df(gaps_df, hgvs_df, annovar_df, snpeff_df, mutalyzer_df, vep_refseq_nomenclature_df, vep_hg19_nomenclature_df)

    # Left join all datasets to the variant list. Write out all rows and all fields  
    jc.write(args.out, gaps_df, hgvs_df, annovar_df, snpeff_df, mutalyzer_df, vep_refseq_nomenclature_df, vep_hg19_nomenclature_df, comparison_df)
    

if __name__ == '__main__':
    main()
