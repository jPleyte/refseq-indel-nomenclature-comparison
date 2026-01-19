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
from rinc.variant_transcript import VariantTranscript
from collections import Counter
from enum import Enum

class NomenclatureTools(Enum):
    HGVS = "hgvs"
    VEP = "vep"
    TFX = "tfx"
    ANNOVAR = "annovar"
    SNPEFF = "snpeff"
    CGD = "cgd"
    VEP_REFSEQ = "vep_refseq"
    VEP_HG19 = "vep_hg19"

class JoinAndCompare(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        
    def _get_all_variant_transcripts(self, dataframes: list[pd.DataFrame]) -> list[VariantTranscript]:
        """
        Return a list of every variant transcript in all of the dataframes 
        """
        join_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        
        merged_df = reduce(lambda left, right: pd.merge(left, right, on=join_cols, how='outer'), dataframes)
        
        variant_transcripts = []
        for row in merged_df[['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']].itertuples(index=False):
            variant_transcripts.append(VariantTranscript(row.chromosome, row.position, row.reference, row.alt, row.cdna_transcript))
            
        return variant_transcripts
    
    def get_nomenclature_df(self, nomenclature_tool: str, csv_file: str):
        """
        Read a file containing variant transcript nomenclature.
        If the file doesn't exist then an empty dataframe is returned.
        """
        df = pd.read_csv(csv_file, dtype=str)
        df.attrs['nomenclature_tool'] = nomenclature_tool
        self._logger.info(f"Read {df.shape[0]} rows for {nomenclature_tool} from {csv_file}")
        return df

    def _get_variant_transcript_row(self, df, chromosome, position, reference, alt, cdna_transcript):
        """
        Query a dataframe for a row with the variant and transcript
        """
        create_mask = lambda df, c, p, r, a, t: (
            (df['chromosome'] == c) & 
            (df['position'] == p) & 
            (df['reference'] == r) & 
            (df['alt'] == a) &
            (df['cdna_transcript'] == t))
        
        df_row = df[create_mask(df, chromosome, position, reference, alt, cdna_transcript)]
        return df_row

        
    def get_comparison_df(self, dataframes: list[pd.DataFrame]):
        """
        Compare dataframes with each other and return a dataframe with comparison score         
        """
        
        # These two lists of dict will be converted to datframes and then merged into one
        pairwise_comparisons = []
        is_of_interest = [] 
        
        all_variants = self._get_all_variant_transcripts(dataframes)
        
        n=0
        total = len(all_variants)
        
        for v in all_variants:
            chromosome = v.chromosome
            position = v.position
            reference = v.reference
            alt = v.alt
            cdna_transcript = v.cdna_transcript
            
            tool_transcript_counter = Counter()
            
            # From each dataource fetch a row match the vairant and transcript  
            rows = []
            for df in dataframes:
                row = self._get_variant_transcript_row(df, chromosome, position, reference, alt, cdna_transcript)
                
                if row.shape[0] > 1:
                    raise ValueError(f"Data source {df.attrs['nomenclature_tool']} has multiple rows for: {chromosome, position, reference, alt, cdna_transcript}")
                if row.empty:
                    self._logger
                    continue
                
                tool_transcript_counter[df.attrs['nomenclature_tool']] += 1
                rows.append(row)
                        
            # Perform pairwise comparison of exon,c.,p. for all datasources that have a matching row                     
            pairwise_comparisons.append(self._get_pairwise_comparison(chromosome, position, reference, alt, cdna_transcript, rows))
            
            # Flag rows that raise an eyebrow 
            is_of_interest.append(self._get_is_of_interest(chromosome, position, reference, alt, cdna_transcript, rows))
            
            n = n + 1
            if n % 10 == 0:
                self._logger.info(f"Processed {n}/{total} transcripts")

        # Create dataframes from the dict objects 
        pairwise_df  = pd.DataFrame(pairwise_comparisons).astype({'pairwise_score': float})
        of_interest_df = pd.DataFrame(is_of_interest)
        
        # Merge the two dataframes    
        join_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        comparison_df = pd.merge(pairwise_df, of_interest_df, on=join_cols, how='left')
            
        assert pairwise_df.shape[0] == of_interest_df.shape[0] == comparison_df.shape[0], "unintended df shape change"
        # if pairwise_df.shape[0] != of_interest_df.shape[0] != comparison_df.shape[0] != "unintended df shape change":
            # raise ValueError("unintended df shape change")

        return comparison_df

    def _get_field_name(self, pattern, row: pd.DataFrame, optional=False):
        """
        Return the name of the field that has c. in it
        """
        field_names = []
        for x in row.columns:
            if pattern in x:
                field_names.append(x)
    
        if len(field_names) > 1: 
            raise ValueError(f"Multiple field names matching {pattern}: {field_names}")
        if not optional and not field_names:
            raise ValueError(f"Unable to find field matching {pattern} in {row.attrs['nomenclature_tool']}")
        elif field_names:
            return field_names[0]
        else:
            return None
    
    def _get_pairwise_comparison(self, chromosome, position, reference, alt, cdna_transcript, rows: list[pd.DataFrame]) -> dict:
        """
        Calculate pairwise comparison score for c., p., and exon from all sources.
        """
        comparison = {'chromosome': chromosome,
                       'position': position,
                       'reference': reference,
                       'alt': alt,
                       'cdna_transcript': cdna_transcript}

        # c. values for comparison
        c_dot_values = []
        
        for row in rows:
            if not row.empty:
                c_dot_values.append(row[self._get_field_name('c_dot', row)].item())
             
        # p. values for comparison        
        p_dot_values = []
        for row in rows:
            if not row.empty:
                p_dot_values.append(row[self._get_field_name('p_dot1', row)].item())    

                     
        # exon vlaues for comparison
        exon_values = []
        for row in rows:
            if not row.empty:
                # Not all dataframes have an exon field
                exon_field_name = self._get_field_name('exon', row, True)
                if exon_field_name:
                    exon_values.append(row[exon_field_name].item())
        
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
            # a score of -3 means only one source provided a value 
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
    
    def _get_is_of_interest(self, chromosome, position, reference, alt, cdna_transcript, rows: list[pd.DataFrame]) -> bool:
        """
        Return a dataframe with a single row "of_interst" that will be true when
        * Any of the g_dot fields don't match
        * Mutalyzer status = Failed
        * Vep's GIVEN_REF and USED_REF don't match
        * Vep's vep.refseq.BAM_EDIT is a value we haven't seen before
        * Vep's vep (either) REFSEQ_MATCH is a value we haven't seen before        
        """
        
        is_of_interest = {'chromosome': chromosome,
                       'position': position,
                       'reference': reference,
                       'alt': alt,
                       'cdna_transcript': cdna_transcript,
                       'of_interest': None 
                       }
        
        interest_codes = []
        
        # Are all of the g. the same? 
        g_dots = set()
        for row in rows:
            g_dot_field = self._get_field_name('g_dot', row, True)
            if g_dot_field:
                g_dots.add(row[g_dot_field].item())
        
        if len(g_dots) > 1:
            self._logger.info(f"Variant g. don't match for {chromosome}-{position}-{reference}-{alt}")
            interest_codes.append('gdot_mismatch') 
                     
        # In the VEP refseq comparison do GIVEN_REF and USED_REF  match?
        vep_refseq_row = self._get_row(NomenclatureTools.VEP_REFSEQ, rows)
        if vep_refseq_row is not None and vep_refseq_row['vep.refseq.GIVEN_REF'].item() != vep_refseq_row['vep.refseq.USED_REF'].item():
            interest_codes.append('vep_ref_difference')
        
        # if the BAM_EDIT is a value we haven't seen before
        vep_refseq_row = self._get_row(NomenclatureTools.VEP_REFSEQ, rows)
        if vep_refseq_row is not None and vep_refseq_row['vep.refseq.BAM_EDIT'].item() not in ['-', 'OK', 'FAILED']:
            interest_codes.append('vep_refseq_unknown_bamedit')
        
        # If REFSEQ_MATCH value we haven't seen before 
        vep_refseq_row = self._get_row(NomenclatureTools.VEP_REFSEQ, rows)
        if vep_refseq_row is not None and vep_refseq_row['vep.refseq.REFSEQ_MATCH'].item() != '-':
            interest_codes.append('vep_refseq_unknown_refseq_match')
        
        # a REFSEQ_MATCH value we haven't seen before
        vep_hg19_row = self._get_row(NomenclatureTools.VEP_HG19, rows)
        if vep_hg19_row is not None and not vep_hg19_row.empty and vep_hg19_row['vep.hg19.REFSEQ_MATCH'].item() != '-':
            interest_codes.append('vep_hg19_refseq_match')
        
        is_of_interest['of_interest'] = '|'.join(interest_codes)
        
        return is_of_interest
        
    def _get_row(self, nomenclature_tool, rows: list):
        """
        Return a row where the  nomenclature_tool attribute matchces the parameter  
        """
        for x in rows:
            if x.attrs['nomenclature_tool'] == nomenclature_tool:
                return x
        
        return None
        
    def write(self, 
              out_file, 
              dataframes: list[pd.DataFrame]):
        """
        Join all the dataframs on the variant+transcript fields and write out a csv with all fields.  
        """        
        join_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        sort_cols = ['pairwise_score', 'chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        sort_orders = [False, True, True, True, True, True]
        
        # Perform outter join on variant+transcript fields         
        merged_df = (
            reduce(lambda left, right: pd.merge(left, right, on=join_cols, how='outer'), dataframes)
            .sort_values(by=sort_cols, ascending=sort_orders)
        )
        
        # There are a lot of fields in the dataframe. Move the important ones up front         
        front_columns = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript', 'pairwise_score']
        
        # Group all the exon columns tobether 
        exons_field_names = list(filter(None, [self._get_field_name('exon', x, True) for x in dataframes]))
        front_columns.extend(exons_field_names)
        
        # Group all the c_dot columns together
        c_dot_field_names = list(filter(None, [self._get_field_name('c_dot', x, True) for x in dataframes]))
        front_columns.extend(c_dot_field_names)
        
        # Group all the p_dot1 columns 
        p_dot1_field_names = list(filter(None, [self._get_field_name('p_dot1', x, True) for x in dataframes]))
        front_columns.extend(p_dot1_field_names)
        
        # Group all the g_dot columns 
        g_dot_field_names = list(filter(None, [self._get_field_name('g_dot', x, True) for x in dataframes]))
        front_columns.extend(g_dot_field_names)
        
        # Group all the pdot3
        p_dot3_field_names = list(filter(None, [self._get_field_name('p_dot3', x, True) for x in dataframes]))
        front_columns.extend(p_dot3_field_names)
        
        # Group the protein transcript columns 
        protein_transcript_field_names = list(filter(None, [self._get_field_name('protein_transcript', x, True) for x in dataframes]))
        front_columns.extend(protein_transcript_field_names)
        
        # Group the gene columns
        gene_field_names = list(filter(None, [self._get_field_name('gene', x, True) for x in dataframes]))
        front_columns.extend(gene_field_names)
        
        # Group genomic region type
        grt_field_names = list(filter(None, [self._get_field_name('genomic_region_type', x, True) for x in dataframes]))
        front_columns.extend(grt_field_names)
        grt_field_names = list(filter(None, [self._get_field_name('grt', x, True) for x in dataframes]))
        front_columns.extend(grt_field_names)
        
        # Group protein variant type
        pvt_field_names = list(filter(None, [self._get_field_name('protein_variant_type', x, True) for x in dataframes]))
        front_columns.extend(pvt_field_names)
        pvt_field_names = list(filter(None, [self._get_field_name('pvt', x, True) for x in dataframes]))
        front_columns.extend(pvt_field_names)
 
        # Gather a list of all columns other than the front_columns
        other_columns = [c for c in merged_df.columns if c not in front_columns]

        merged_df.to_csv(out_file, index=False, encoding='utf-8', columns=front_columns+other_columns)
        
        self._logger.info(f"Wrote data frame with {merged_df.shape[0]} rows and {merged_df.shape[1]} columns  to {out_file}")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Read Annovar multianno file, extract values and write to new csv')

    parser.add_argument('--variants', help='Variants (csv)', required=True)
    
    parser.add_argument('--hgvs_nomenclature', help='hgvs/uta values (csv)', required=False)
    parser.add_argument('--tfx_nomenclature', help="Optional Transcript Effect/tfx values (csv)", required=False)
    parser.add_argument('--cgd_nomenclature', help='CGD values (csv)', required=False)
    
    parser.add_argument('--annovar_nomenclature', help='Annovar values (csv)', required=True)    
    parser.add_argument('--snpeff_nomenclature', help='SnpEff values (csv)', required=True)    
    parser.add_argument('--vep_refseq_nomenclautre', help='mutalyzer values (csv)', required=True)
    parser.add_argument('--vep_hg19_nomenclature', help='mutalyzer values (csv)', required=True)
    
    parser.add_argument("--out", help="output file (csv)", required=True)

    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    jc = JoinAndCompare()
    
    dataframes = []
    
    # The hgvs and tfx dataframes are optional
    if args.hgvs_nomenclature:        
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.HGVS, args.hgvs_nomenclature))
    if args.tfx_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.TFX, args.tfx_nomenclature))
    if args.cgd_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.CGD, args.cgd_nomenclature))
        
    # In the future i will make this optional but for now they're always being generated 
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.ANNOVAR, args.annovar_nomenclature))    
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.SNPEFF, args.snpeff_nomenclature))        
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.VEP_REFSEQ, args.vep_refseq_nomenclautre))    
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.VEP_HG19, args.vep_hg19_nomenclature))
    
    # Compare exon, c., and p. from all datasources. 
    comparison_df = jc.get_comparison_df(dataframes)
    dataframes.append(comparison_df)
    
    # Left join all datasets to the variant list. Write out all rows and all fields
    jc.write(args.out, dataframes)
    

if __name__ == '__main__':
    main()
