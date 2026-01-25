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
import itertools
import numpy as np
from openpyxl.styles import Font

class NomenclatureTools(Enum):
    HGVS = "hgvs"
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
    
    def get_nomenclature_df(self, nomenclature_tool: str, csv_file: str, remove_label: str):
        """
        Read a file containing variant transcript nomenclature.
        If the file doesn't exist then an empty dataframe is returned.
        
        The dataframe columns all have identifiers in the (eg tfx.exon or p_dot1.annovar) but they get in the way 
        here so the `removelabel` string is removed from the columns. 
        """
        df = pd.read_csv(csv_file, dtype=str)
        
        index_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        df.set_index(index_cols, inplace=True, drop=False)

        df.columns = [x.replace(remove_label, '') for x in df.columns]
        df.columns = pd.MultiIndex.from_product([[nomenclature_tool.value], df.columns])
        
        df.attrs['nomenclature_tool'] = nomenclature_tool.value
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

    def _get_merged_df(self, dataframes: list[pd.DataFrame]):
        """
        Merge all the dataframes into one. 
        I don't want to inner join them all because that is very limiting.
        I don't want to outter join them all because that will give me too much.  
        So i inner join on certain high qualty sources like annovar, tfx, and vep
        And then left join the others in.  
               
        If i wanted to inner join them all i would use: pd.concat(dataframes, axis=1, join='inner')
          
        """
        inner_join_dfs = []  
        outer_join_dfs = []
        
        for x in dataframes:
            if x.attrs['nomenclature_tool'] == NomenclatureTools.HGVS.value:
                outer_join_dfs.append(x)
            elif x.attrs['nomenclature_tool'] == NomenclatureTools.ANNOVAR.value:
                inner_join_dfs.append(x)
            elif x.attrs['nomenclature_tool'] == NomenclatureTools.CGD.value:
                outer_join_dfs.append(x)
            elif x.attrs['nomenclature_tool'] == NomenclatureTools.SNPEFF.value:
                outer_join_dfs.append(x)
            elif x.attrs['nomenclature_tool'] == NomenclatureTools.TFX.value:
                inner_join_dfs.append(x)
            elif x.attrs['nomenclature_tool'] == NomenclatureTools.VEP_HG19.value:
                inner_join_dfs.append(x)
            elif x.attrs['nomenclature_tool'] == NomenclatureTools.VEP_REFSEQ.value:
                inner_join_dfs.append(x)
            else:
                raise ValueError(f"Unknown tool: {x.attrs['nomenclature_tool']}")
                
        inner_df = pd.concat(inner_join_dfs, axis=1, join='inner')
        
        aux_df = pd.concat(outer_join_dfs, axis=1, join='outer')
                
        merged_df = inner_df.join(aux_df, how='left')
        return merged_df    
        
    def get_comparison_df(self, dataframes: list[pd.DataFrame]):
        """
        Compare dataframes with each other and return a dataframe with comparison score         
        """
        merged_df = self._get_merged_df(dataframes)
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['g_dot'], new_field_name='g_dot_concordance')
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['c_dot'], new_field_name='c_dot_concordance')
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['p_dot1'], new_field_name='p_dot_cordance')
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['c_dot', 'p_dot1'], new_field_name='c+p_concordance')
        merged_df = self._calculate_pairwise_score(merged_df, dataframes, fields_to_compare=['exon', 'c_dot', 'p_dot1'], new_field_name='c+p+exon_concordance')
        merged_df = self._add_vep_refseq_ref_mismatch_field(merged_df, new_field_name='vep_refseq_mismatch')
        
        # Add a field indicating when cgd&annovar agree but disagree with tfx&vepHg19 on c_dot
        merged_df = self._add_consensus_conflict_field(merged_df, ['cgd', 'annovar'], ['tfx', 'vep_hg19'], ['c_dot'], 'ca_vs_tvv_conflict')
        
                        
        return merged_df

    def _calculate_pairwise_score(self, merged_df: pd.DataFrame, tool_dataframes: list[pd.DataFrame], fields_to_compare: list[str], new_field_name):
        """
        Compare the fields_to_compare   
        """
        tool_dataframe_names = [x.attrs['nomenclature_tool'] for x in tool_dataframes]
        
        pairs = list(itertools.combinations(tool_dataframe_names, 2))
        
        total_matches = pd.Series(0, index=merged_df.index)
                
        possible_comparisons = pd.Series(0, index=merged_df.index)
        
        for t1, t2 in pairs:
            
            t1_has_fields = all((t1, f) in merged_df.columns for f in fields_to_compare)
            t2_has_fields = all((t2, f) in merged_df.columns for f in fields_to_compare)
            if not (t1_has_fields and t2_has_fields):
                self._logger.info(f"Unable to calculate pairwise score of {fields_to_compare} because {t1} or {t2} does not have that one of the fields.")
                continue
                 
            # 1. Determine if BOTH tools have data for ALL fields in the list
            # We start with True and 'AND' it with each field's presence
            both_have_all_data = pd.Series(True, index=merged_df.index)
            for field in fields_to_compare:
                both_have_all_data &= (merged_df[t1, field].notna() & merged_df[t2, field].notna())
            
            # 2. Determine if ALL fields match exactly between the two tools
            # We start with True and 'AND' it with each field's comparison
            all_fields_match = pd.Series(True, index=merged_df.index)
            for field in fields_to_compare:
                all_fields_match &= (merged_df[t1, field] == merged_df[t2, field])
            
            # 3. A 'Point' is awarded only if they match AND both had data
            match_mask = all_fields_match & both_have_all_data
            
            # 3. Add to counters
            total_matches += match_mask.astype(int)
            possible_comparisons += both_have_all_data.astype(int)
            
        # 4. Calculate Score: matches / comparisons 
        # (avoid division by zero with .replace)        
        merged_df['scores', new_field_name] = total_matches / possible_comparisons.replace(0, np.nan)
        return merged_df
    
        
    def _get_field_name(self, pattern, df: pd.DataFrame, optional=False):
        """
        Return the name of the field that has c. in it
        """
        field_names = []
        for x in df.columns:
            if pattern in x:
                field_names.append(x)
    
        if len(field_names) > 1: 
            raise ValueError(f"Multiple field names matching {pattern}: {field_names}")
        if not optional and not field_names:
            raise ValueError(f"Unable to find field matching {pattern} in {df.attrs['nomenclature_tool']}")
        elif field_names:
            return field_names[0]
        else:
            return None
    
    def _add_vep_refseq_ref_mismatch_field(self, merged_df: pd.DataFrame, new_field_name: str):
        """
        Add a field that is 1 when vepRefSeq GIVEN_REF != USED_REF
        """
        mismatch_mask = (
            (merged_df[(NomenclatureTools.VEP_REFSEQ.value, 'GIVEN_REF')] != merged_df[(NomenclatureTools.VEP_REFSEQ.value, 'USED_REF')]) & 
            merged_df[(NomenclatureTools.VEP_REFSEQ.value, 'GIVEN_REF')].notna() & 
            merged_df[(NomenclatureTools.VEP_REFSEQ.value, 'USED_REF')].notna()
            )
        
        merged_df[('scores', new_field_name)] = mismatch_mask.astype(int)
        return merged_df
        
    def _get_row(self, nomenclature_tool, rows: list):
        """
        Return a row where the  nomenclature_tool attribute matchces the parameter  
        """
        for x in rows:
            if x.attrs['nomenclature_tool'] == nomenclature_tool:
                return x
        
        return None
        
    def _write_sheet_raw(self, writer, workbook, df, sheet_name):
        """
        Add the raw data worksheet 
        """
        df.to_excel(writer, sheet_name=sheet_name)
        worksheet = writer.sheets[sheet_name]
        
        bold_fmt = workbook.add_format({'bold': True, 'bg_color': '#D3D3D3', 'border': 1})
        worksheet.set_row(0, None, bold_fmt)
        worksheet.set_row(1, None, bold_fmt)
        
        worksheet.freeze_panes(2, 0)
        worksheet.set_column('A:ZZ', 18)

    def _write_sheet_comparison(self, writer, workbook, flat_df: pd.DataFrame, sheet_name):
        """
        Format the data to make it easier to read  
        """
        flat_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Bold the single header row in the new sheet
        summary_sheet = writer.sheets[sheet_name]
        
        bold_fmt = workbook.add_format({'bold': True, 'bg_color': '#D3D3D3', 'border': 1})
        summary_sheet.set_row(0, None, bold_fmt)
        
        summary_sheet.freeze_panes(1, 0)
        summary_sheet.set_column('A:E', 12) # Genomic coords
        summary_sheet.set_column('F:ZZ', 25) # Nomenclature fields
        
    def _write_sheet_summary(self, writer, workbook, flat_df, sheet_name):
        """
        """
        # 1. Identify all your score columns
        score_cols = [c for c in flat_df.columns if 'scores_' in c]
        
        summary_data = []
        
        for col in score_cols:
            # Get the data for this specific score
            data = flat_df[col]
            
            # How many variants had enough data to actually produce a score?
            total_comparable = data.notna().sum()
            
            # How many of those were perfect matches (score == 1.0)?
            perfect_matches = (data == 1.0).sum()
            
            # Calculate percentage
            percentage = (perfect_matches / total_comparable * 100) if total_comparable > 0 else 0
            
            summary_data.append({
                'Metric': col.replace('scores_', ''),
                'Perfect Concordance (n)': perfect_matches,
                'Total Comparable Variants': total_comparable,
                'Concordance Rate (%)': round(percentage, 2)
            })
        
        # 2. Create the Summary DataFrame
        metrics_df = pd.DataFrame(summary_data)
        
        # 3. Write to a new sheet in your Excel file
        # (Add this inside your 'with pd.ExcelWriter...' block)
        metrics_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # 4. Optional: Add a little formatting to the Metrics sheet
        ws_metrics = writer.sheets[sheet_name]
        ws_metrics.set_column('A:A', 35) # Make the Metric name wide
        ws_metrics.set_column('B:D', 20) # Center the numbers
        
    def _get_flattened_data_frame(self, df):
        """
        """
        # We reset the index to turn the 5 genomic levels (Chr, Pos, etc.) into regular columns        
        flat_df = df.copy().reset_index()
        
        # 2. Force everything to string, handling MultiIndex tuples specifically
        new_cols = []
        for c in flat_df.columns:
            if isinstance(c, tuple):
                # Join the tuple parts with underscore, filter out empty strings
                clean_col = "_".join([str(part) for part in c if str(part).strip()])
                new_cols.append(clean_col)
            else:
                new_cols.append(str(c))
        
        
        flat_df.columns = new_cols
        
        # Identify columns index columns
        index_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
    
        # every dataframe brought in its own index column (eg cgd_chromosome) and we don't need to see those repeated.  
        redundant_suffixes = [
            'chromosome', 'position', 'reference', 'alt', 'cdna_transcript'
        ] 
           
        # nom_cols are all the not index "nomenclature" and other fields 
        nom_cols = [c for c in flat_df.columns 
                    if c not in index_cols
                    and not any(c.endswith(suffix) for suffix in redundant_suffixes)
                    ]
                                
        field_order = ['c_dot', 'p_dot1', 'exon', 'g_dot']

        # 4. The Sort
        sorted_nom_cols = sorted(
            nom_cols, 
            key=lambda x: (
                next((i for i, f in enumerate(field_order) if x.endswith(f)), len(field_order)),
                x
            )
        )
        
        final_column_order = index_cols + sorted_nom_cols
        return flat_df[final_column_order]        
        
    def _add_consensus_conflict_field(self, merged_df: pd.DataFrame, list_a: list[str], list_b: list[str], fields: list[str], new_field_name):
        """
        Identifies rows where ListA tools agree, ListB tools agree, but Group A disagrees with Group B.
        """
        # Initialize a mask that starts as True for all rows
        # We will 'AND' this with our conditions
        final_mask = pd.Series(True, index=merged_df.index)
    
        for field in fields:
            # --- 1. Check Internal Consensus for List A ---
            # All must be non-null and equal to the first tool in the list
            a_first = list_a[0]
            a_consensus = pd.Series(True, index=merged_df.index)
            
            for tool in list_a:
                # Must have data
                a_consensus &= merged_df[tool, field].notna()
                # Must match the first tool
                a_consensus &= (merged_df[tool, field] == merged_df[a_first, field])
            
            # --- 2. Check Internal Consensus for List B ---
            b_first = list_b[0]
            b_consensus = pd.Series(True, index=merged_df.index)
            
            for tool in list_b:
                b_consensus &= merged_df[tool, field].notna()
                b_consensus &= (merged_df[tool, field] == merged_df[b_first, field])
                
            # --- 3. Check for Disagreement between Group A and Group B ---
            groups_disagree = (merged_df[a_first, field] != merged_df[b_first, field])
            
            # Combine for this specific field
            field_criteria = a_consensus & b_consensus & groups_disagree
            
            # Update the final mask (Row must meet criteria for ALL fields, e.g., c_dot AND p_dot)
            final_mask &= field_criteria
    
        self._logger.info(f"Consensus conflict for {list_a} vs {list_b} has {final_mask.sum()} rows")
        
        # Add a boolean field indicating rows that match the criteria
        merged_df[('scores', new_field_name)] = final_mask
        return merged_df
        

    def write(self, 
              out_file, 
              comparison_df: pd.DataFrame):
        """
        Write the dataframe with its new comparison fields to file
        #merged_df.columns = [f"{source}_{field}" for source, field in merged_df.columns]
        #merged_df.to_csv("/tmp/merged_pairwise_df.csv", index=True)  
        """                        
        #sort_cols = ['c+p+exon_concordance', 'chromosome', 'position', 'reference', 'alt']
        #sort_orders = [False, True, True, True, True]
        comparison_df.sort_values(
            by=('scores', 'c+p+exon_concordance'), 
            ascending=False, 
            inplace=True
            )
        
        self._logger.info("Creating workbook")
        
        with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
            workbook  = writer.book
            
            flatened_df = self._get_flattened_data_frame(comparison_df)
            
            self._write_sheet_summary(writer, workbook, flatened_df, 'Summary')            
            self._write_sheet_comparison(writer, workbook, flatened_df, 'Comparison')
            self._write_sheet_raw(writer, workbook, comparison_df, 'Raw')
            

        self._logger.info(f"Wrote data frame with {comparison_df.shape[0]} rows and {comparison_df.shape[1]} columns to {out_file}")
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Read Annovar multianno file, extract values and write to new csv')

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
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.HGVS, args.hgvs_nomenclature, 'hu.'))
    if args.tfx_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.TFX, args.tfx_nomenclature, '.tfx'))
    if args.cgd_nomenclature:
        dataframes.append(jc.get_nomenclature_df(NomenclatureTools.CGD, args.cgd_nomenclature, '.cgd'))
        
    # In the future i will make this optional but for now they're always being generated 
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.ANNOVAR, args.annovar_nomenclature, '.annovar'))    
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.SNPEFF, args.snpeff_nomenclature, 'snpeff.'))
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.VEP_REFSEQ, args.vep_refseq_nomenclautre, 'vep.refseq.'))    
    dataframes.append(jc.get_nomenclature_df(NomenclatureTools.VEP_HG19, args.vep_hg19_nomenclature, 'vep.hg19.'))
    
    # Compare exon, c., and p. from all datasources. 
    comparison_df = jc.get_comparison_df(dataframes)
    
    # Left join all datasets to the variant list. Write out all rows and all fields
    jc.write(args.out, comparison_df)
    

if __name__ == '__main__':
    main()
