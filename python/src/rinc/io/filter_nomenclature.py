'''
Created on Jan 21, 2026

@author: pleyte
'''
import pandas as pd
import argparse
import sys
import os

def main():    
    parser = argparse.ArgumentParser(description="Filter nomenclature CSV based on a whitelist CSV.")
    parser.add_argument('--input', required=True, help="The tool output CSV to be filtered.")
    parser.add_argument('--filter', required=True, help="The filter/whitelist CSV.")
    parser.add_argument('--output', required=True, help="Path to save the filtered CSV.")
    
    args = parser.parse_args()
    
    # 1. Load the input data
    if not os.path.exists(args.input):
        sys.exit(f"Error: Input file {args.input} not found.")
    df_input = pd.read_csv(args.input, dtype=str)

    # 2. Load the filter data
    if not os.path.exists(args.filter):
        sys.exit(f"Error: Filter file {args.filter} not found.")
    df_filter = pd.read_csv(args.filter, dtype=str)

    # 3. Define the columns to match on
    # These must exist in both files
    match_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']

    # 4. Perform the filter using an inner merge
    # This keeps only rows where the combination of these 5 columns 
    # exists in both the input and the filter file.
    filtered_df = pd.merge(
        df_input, 
        df_filter[match_cols], 
        on=match_cols, 
        how='inner'
    )

    # 5. Write the output
    filtered_df.to_csv(args.output, index=False)
    
    print(f"Filtering complete. Kept {len(filtered_df)} out of {len(df_input)} rows.")    

if __name__ == "__main__":
    main()