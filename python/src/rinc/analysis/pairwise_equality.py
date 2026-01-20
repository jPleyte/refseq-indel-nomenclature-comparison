'''
Created on Jan 18, 2026

@author: pleyte
'''
from rinc.util.log_config import LogConfig
import argparse
import logging.config
import pandas as pd
from itertools import combinations
    
class PairwiseEqualityAnalysis(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._data_frames = []
        self._analyses: list[dict] = []
    
    def add_datasource(self, name, path):
        """
        """
        df = pd.read_csv(path, dtype=str)
        df.attrs['nomenclature_tool'] = name
        
        self._logger.info(f"Read {df.shape[0]} rows for tool {name} from {path}")
        self._data_frames.append(df)
    
    def analyze(self):
        """
        """
        for (df_a, df_b) in list(combinations(self._data_frames, 2)):
            analysis = self._analyze(df_a, df_b)
            self._analyses.append(analysis)
            self._logger.info(f"Completed {analysis['row_count']} row analysis of {analysis['source_a']} and {analysis['source_b']}")
            
            
    def _analyze(self, df_a: pd.DataFrame, df_b: pd.DataFrame):
        analysis = {}
        analysis['source_a'] = df_a.attrs['nomenclature_tool']
        analysis['source_b'] = df_b.attrs['nomenclature_tool']
        
        columns = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        combined_df = pd.merge(df_a, df_b, on=columns, how='inner')
        
        analysis['row_count'] = combined_df.shape[0]
        
        analysis['g_dot_equal'] = self._get_count_equal_field('g_dot', df_a, df_b, combined_df, optional=True)
        analysis['c_dot_equal'] = self._get_count_equal_field('c_dot', df_a, df_b, combined_df)
        analysis['p_dot_equal'] = self._get_count_equal_field('p_dot1', df_a, df_b, combined_df)
        
        analysis['c+p_dot_equal'] = self._get_count_c_and_p_equal(None, df_a, df_b, combined_df)
        analysis['cp_equal/refseq'] = self._get_count_c_and_p_equal('NM', df_a, df_b, combined_df)
        analysis['cp_equal/ccds'] = self._get_count_c_and_p_equal('CCDS', df_a, df_b, combined_df)
        
        return analysis
        
    def _get_count_c_and_p_equal(self, transcript_prefix:str, df_a: pd.DataFrame, df_b: pd.DataFrame, combined_df: pd.DataFrame) -> str:
        """
        Count the number of rows where two datasources have c. and p. in common specific to RefSeq or CCDS transcripts
        """
        df_a_c_dot_field = self._get_field_name_matching('c_dot', df_a)
        df_a_p_dot_field = self._get_field_name_matching('p_dot1', df_a)
        cols_a = [df_a_c_dot_field, df_a_p_dot_field]
        
        df_b_c_dot_field = self._get_field_name_matching('c_dot', df_b)
        df_b_p_dot_field = self._get_field_name_matching('p_dot1', df_b)
        cols_b = [df_b_c_dot_field, df_b_p_dot_field]

        # Create a mask for c. and p. equality
        nomenclature_match = (combined_df[cols_a].values == combined_df[cols_b].values).all(axis=1)
        
        # Create a mask for transcripts that start with the transcript_prefix
        if transcript_prefix:
            transcript_filter = combined_df['cdna_transcript'].str.startswith(transcript_prefix, na=False)
            cnt_with_prefix = transcript_filter.sum()
            cnt_with_prefix_and_cp_match = (nomenclature_match & transcript_filter).sum()            
        else:
            cnt_with_prefix = combined_df.shape[0]
            cnt_with_prefix_and_cp_match = (nomenclature_match).sum()
            return f"{cnt_with_prefix_and_cp_match} ({round(((cnt_with_prefix_and_cp_match / cnt_with_prefix)*100),1)}%)"
                     
        if cnt_with_prefix == 0:
            return None
        
        
        assert cnt_with_prefix_and_cp_match <= cnt_with_prefix 
        
        return f"{round(((cnt_with_prefix_and_cp_match / cnt_with_prefix)*100),1)}%"


    def _get_count_equal_field(self, field:str, df_a: pd.DataFrame, df_b: pd.DataFrame, combined_df: pd.DataFrame, optional=False) -> int:
        """
        Count the number of rows where two datasources have a single field (like c.) in common 
        """
        df_a_field = self._get_field_name_matching(field, df_a, optional)
        df_b_field = self._get_field_name_matching(field, df_b, optional)
        
        if optional and (not df_a_field or not df_b_field):
            return None
        
        return (combined_df[df_a_field] == combined_df[df_b_field]).sum()
    
    def _get_field_name_matching(self, field: str, df: pd.DataFrame, optional=False):
        """
        Find a field name in the dataframe that partially matches the field parameter (eg 'vep.c_dot' for 'c_dot') 
        """
        matches = [item for item in df.columns if field in item]
        
        if len(matches) == 0 and optional:
            return None
        elif len(matches) == 0:
            raise ValueError(f"No fields matching {field} in field list for {df.attrs['nomenclature_tool']}")
        elif len(matches) > 1:
            raise ValueError(f"More than one field matches {field} in field list for {df.attrs['nomenclature_tool']}: {matches}")
            
        return matches[0]
        
    def write(self, output_filename: str):
        """
        Write analysis to file
        """
        headers_source = ['source_a', 'source_b', 'row_count', 'g_dot_equal', 'c_dot_equal', 'p_dot_equal', 'c+p_dot_equal', 'cp_equal/refseq', 'cp_equal/ccds']
        
        dtype_schema = {
            'source_a': 'str', 
            'source_b': 'str',
            'row_count': 'Int64',
            'g_dot_equal': 'Int64',
            'c_dot_equal': 'Int64',
            'p_dot_equal': 'Int64',
            'c+p_dot_equal': 'str',
            'cp_equal/refseq': 'str',
            'cp_equal/ccds': 'str'
        }
        analyses_df = pd.DataFrame(self._analyses).astype(dtype_schema)
        
        # add a new field that is the sum of the *_dot_equal fields
        analyses_df['x_dot_equal_match_rate'] = (
            analyses_df[['c_dot_equal', 'p_dot_equal']].sum(axis=1) / analyses_df['row_count']
            )
        # analyses_df['x_dot_equal_match_sum'] = analyses_df[['c_dot_equal', 'p_dot_equal']].sum(axis=1)  / analyses_df.sum()
        
        analyses_df.sort_values(by='x_dot_equal_match_rate', ascending=False) \
            .to_csv(output_filename, index=False, encoding='utf-8', columns=headers_source, na_rep='')
                  
        self._logger.info("Wrote pairwise equality analysis to "+output_filename)


def _parse_args():
    parser = argparse.ArgumentParser(description='Read Annovar multianno file, extract values and write to new csv')

    parser.add_argument('--nomenclature', action='append', nargs=2, metavar=('NAME', 'PATH'), help="Pass nickname followed by path to CSV", required=True)    
    parser.add_argument("--out", help="output file (csv)", dest="output", required=True)

    parser.add_argument("--version", action="version", version="0.0.1")

    return parser.parse_args()

    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    pea = PairwiseEqualityAnalysis()
    
    if args.nomenclature:
        for name, path in args.nomenclature:
            pea.add_datasource(name, path)
    
    pea.analyze()
    pea.write(args.output)
    

if __name__ == '__main__':
    main()