'''
Created on Jan 27, 2026

@author: pleyte
'''
import argparse
import logging.config
from rinc.util.log_config import LogConfig
import pandas as pd

def _get_join_dataframes(gff_gaps_file, uta_gaps_file, out_inner_join_csv=None):
    """
    Read the list of accessions pulled from the gff  and the list of accessions from UTA
    Join them on their refseq accession and return the resulting df
    """
    gff_df = pd.read_parquet(gff_gaps_file)
    print(f"Read {gff_df.shape[0]} rows from {gff_gaps_file}")
    
    gff_df = gff_df[gff_df['gap'].notna() & (gff_df['gap'] != '')]
    print(f"Reduced gff df to {gff_df.shape[0]} rows having gap info")
    
    uta_df = pd.read_csv(uta_gaps_file)
    print(f"Read {uta_df.shape[0]} rows from {uta_gaps_file}")
    
    uta_df = uta_df[uta_df['cigar'].notna() & (uta_df['cigar'] != '')]
    print(f"Reduced uta df to {uta_df.shape[0]} rows having gap info")
    
    if out_inner_join_csv:
        _write_out_inner_join_csv(out_inner_join_csv, gff_df, uta_df)
        
    combined_df = pd.merge(
        gff_df, 
        uta_df, 
        on='accession', 
        how='outer', 
        suffixes=('_gff', '_uta')
    )
    
    print(f"Outter joined the two dataframes to produce {combined_df.shape[0]} rows")
    
    # Rename fields and only return the fields we want
    return (combined_df.rename(columns={
        'gap': 'gff_cigars',
        'cigar': 'uta_cigars'
        })
        [['accession', 'gff_cigars', 'uta_cigars']]
    )

def _write_out_inner_join_csv(out_inner_join_csv, gff_df, uta_df):
    """
    For debugging, write out a csv that shows us which transcript gaps the gff and uta have in common 
    by using an inner join to join the two dfs.   
    """
    combined_df = pd.merge(
        gff_df, 
        uta_df, 
        on='accession', 
        how='inner', 
        suffixes=('_gff', '_uta')
    )
    
    # Rename fields and only return the fields we want
    final_df = (combined_df.rename(columns={
        'gap': 'gff_cigars',
        'cigar': 'uta_cigars'
        })
        [['accession', 'gff_cigars', 'uta_cigars']]
    )
    
    final_df.to_csv(out_inner_join_csv)
    print(f"Wrote {final_df.shape[0]} rows resulting from inner join of gff and uta dfs to {out_inner_join_csv}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Query UTA for alignment differences')
    parser.add_argument("--version", action="version", version="0.0.1")    
    parser.add_argument("--gff_gaps", help="Transcripts with gaps identified using gff", required=True)
    parser.add_argument("--uta_gaps", help="Transcripts with gaps identified using UTA db", required=True)
    parser.add_argument("--out_inner_join_csv", help="Optional csv resulting from inner join, for additional info", required=False)
    parser.add_argument("--out_csv", help="output file (csv)", required=True)
    args = parser.parse_args()
    return args
    
def main():
    logging.config.dictConfig(LogConfig().stdout_config)
    
    args = _parse_args()
    
    combined_df = _get_join_dataframes(args.gff_gaps, args.uta_gaps, args.out_inner_join_csv)
    combined_df.to_csv(args.out_csv, index=False)
    print(f"Wrote {combined_df.shape[0]} rows to {args.out_csv}")

    
if __name__ == '__main__':
    main()
