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
        df.to_excel(writer, sheet_name='Raw')
        worksheet = workbook['Raw']
        bold_font = Font(bold=True)
        for row in worksheet.iter_rows(min_row=1, max_row=2):
            for cell in row:
                cell.font = bold_font
        
        worksheet.freeze_panes = 'A3'

    def _write_sheet_comparison(self, writer, workbook, df, sheet_name):
        """
        Format the data to make it easier to read  
        """
        human_fields = ['c_dot', 'p_dot1', 'exon', 'g_dot']
        selected_columns = [
            col for col in df.columns if col[1] in human_fields and col[0] != 'scores'
        ]
        
        # We reset the index to turn the 5 genomic levels (Chr, Pos, etc.) into regular columns
        summary_df = df[selected_columns].copy()
        summary_df = summary_df.reset_index()
        
        # 4. Flatten the headers from ('tool', 'field') to 'tool_field'
        # For the index columns (like 'chromosome'), they won't have a second level, 
        # so we handle them gracefully.
        new_headers = []
        for col in summary_df.columns:
            if isinstance(col, tuple) and col[1] != '':
                new_headers.append(f"{col[0]}_{col[1]}")
            else:
                # This handles the index columns like 'chromosome', 'position', etc.
                new_headers.append(col[0] if isinstance(col, tuple) else col)
        
        summary_df.columns = new_headers
        summary_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Bold the single header row in the new sheet
        summary_sheet = writer.book[sheet_name]
        for cell in summary_sheet[1]: # Row 1
            cell.font = Font(bold=True)
        
        summary_sheet.freeze_panes = 'A2' # Freeze just the header
        
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
        
        # with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
        with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
            workbook  = writer.book
            
            self._write_sheet_raw(writer, workbook, comparison_df, 'Raw')
            self._write_sheet_comparison(writer, workbook, comparison_df, 'Comparison')

        self._logger.info("jDebug3")
        #comparison_df.to_excel(writer, sheet_name='Comparison')
        #self._logger.info("jDebug4")

        # join_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        # sort_cols = ['pairwise_score', 'chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        #
        #
        #
        # # There are a lot of fields in the dataframe. Move the important ones up front         
        # front_columns = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript', 'pairwise_score']
        #
        # # Group all the exon columns tobether 
        # exons_field_names = list(filter(None, [self._get_field_name('exon', x, True) for x in dataframes]))
        # front_columns.extend(exons_field_names)
        #
        # # Group all the c_dot columns together
        # c_dot_field_names = list(filter(None, [self._get_field_name('c_dot', x, True) for x in dataframes]))
        # front_columns.extend(c_dot_field_names)
        #
        # # Group all the p_dot1 columns 
        # p_dot1_field_names = list(filter(None, [self._get_field_name('p_dot1', x, True) for x in dataframes]))
        # front_columns.extend(p_dot1_field_names)
        #
        # # Group all the g_dot columns 
        # g_dot_field_names = list(filter(None, [self._get_field_name('g_dot', x, True) for x in dataframes]))
        # front_columns.extend(g_dot_field_names)
        #
        # # Group all the pdot3
        # p_dot3_field_names = list(filter(None, [self._get_field_name('p_dot3', x, True) for x in dataframes]))
        # front_columns.extend(p_dot3_field_names)
        #
        # # Group the protein transcript columns 
        # protein_transcript_field_names = list(filter(None, [self._get_field_name('protein_transcript', x, True) for x in dataframes]))
        # front_columns.extend(protein_transcript_field_names)
        #
        # # Group the gene columns
        # gene_field_names = list(filter(None, [self._get_field_name('gene', x, True) for x in dataframes]))
        # front_columns.extend(gene_field_names)
        #
        # # Group genomic region type
        # grt_field_names = list(filter(None, [self._get_field_name('genomic_region_type', x, True) for x in dataframes]))
        # front_columns.extend(grt_field_names)
        # grt_field_names = list(filter(None, [self._get_field_name('grt', x, True) for x in dataframes]))
        # front_columns.extend(grt_field_names)
        #
        # # Group protein variant type
        # pvt_field_names = list(filter(None, [self._get_field_name('protein_variant_type', x, True) for x in dataframes]))
        # front_columns.extend(pvt_field_names)
        # pvt_field_names = list(filter(None, [self._get_field_name('pvt', x, True) for x in dataframes]))
        # front_columns.extend(pvt_field_names)
        #
        # # Gather a list of all columns other than the front_columns
        # other_columns = [c for c in merged_df.columns if c not in front_columns]
        #
        # merged_df.to_csv(out_file, index=False, encoding='utf-8', columns=front_columns+other_columns)
        
        self._logger.info(f"Wrote data frame with {comparison_df.shape[0]} rows and {comparison_df.shape[1]} columns to {out_file}")
        
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
